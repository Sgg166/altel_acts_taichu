#include "TelFW.hh"
#include "myrapidjson.h"

#include "gl.h"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "TelEvent.hpp"

class glfw_test{
public:
  std::string geometry_path;
  JsonDocument jsd_trans;
  std::unique_ptr<altel::TelGL> telgl;
  glm::vec3 cameraPos;
  glm::vec3 worldCenter;
  glm::vec3 cameraUp;
  float deltaTime = 0.0f; // time between current frame and last frame
  float lastFrame = 0.0f;

  glfw_test(const std::string& path_geometry)
    :geometry_path(path_geometry), telgl(nullptr),
     m_size_ring(10), m_vec_ring_ev(m_size_ring)
  {

  }

  int beginHook(GLFWwindow* window){
    std::string jsstr_geo = altel::TelGL::readFile(geometry_path);
    JsonDocument jsd_geo;
    jsd_geo = altel::TelGL::createJsonDocument(jsstr_geo);
    jsd_trans = altel::TelGL::createTransformExample();

    telgl.reset(new altel::TelGL(JsonValue()));
    telgl->updateGeometry(jsd_geo);

    cameraPos   = glm::vec3(0.0f, 0.0f,  -1000.0f);
    worldCenter = glm::vec3(0.0f, 0.0f,  0.0f);
    cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);
    return 0;
  }
  int clearHook(GLFWwindow* window){
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
    int width;
    int height;
    glfwGetFramebufferSize(window, &width, &height);
    float currentWinRatio = width / (float) height;
    // if(glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
    //   return 0;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
      glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
      cameraPos -= (cameraPos - worldCenter) * (deltaTime * 1.0f);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
      cameraPos += (cameraPos - worldCenter) * (deltaTime * 1.0f);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
      cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
      cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
      cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
      cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);

    jsd_trans["trans"]["lookat"]["eye"]["x"]= cameraPos[0];
    jsd_trans["trans"]["lookat"]["eye"]["y"]= cameraPos[1];
    jsd_trans["trans"]["lookat"]["eye"]["z"]= cameraPos[2];
    jsd_trans["trans"]["lookat"]["center"]["x"]= worldCenter[0];
    jsd_trans["trans"]["lookat"]["center"]["y"]= worldCenter[1];
    jsd_trans["trans"]["lookat"]["center"]["z"]= worldCenter[2];
    jsd_trans["trans"]["lookat"]["up"]["x"]= cameraUp[0];
    jsd_trans["trans"]["lookat"]["up"]["y"]= cameraUp[1];
    jsd_trans["trans"]["lookat"]["up"]["z"]= cameraUp[2];
    jsd_trans["trans"]["persp"]["fov"]= 60;
    jsd_trans["trans"]["persp"]["ratio"]= currentWinRatio;
    jsd_trans["trans"]["persp"]["near"]= 0.1;
    jsd_trans["trans"]["persp"]["far"]= 2000;
    telgl->updateTransform(jsd_trans);
    return 1;
  }

  uint64_t m_size_ring;
  std::vector<std::shared_ptr<altel::TelEvent>> m_vec_ring_ev;
  std::shared_ptr<altel::TelEvent> m_ring_end;//nullptr
  std::atomic_uint64_t m_count_ring_write{1};
  std::atomic_uint64_t m_count_ring_read{1};
  std::atomic_uint64_t m_hot_p_read{0};

  void clearBufferEvent(){ // by write thread
    uint64_t count_ring_read = m_count_ring_read;
    m_count_ring_write = count_ring_read;
  }

  void pushBufferEvent(std::shared_ptr<altel::TelEvent> df){ //by write thread
    uint64_t next_p_ring_write = m_count_ring_write % m_size_ring;
    if(next_p_ring_write == m_hot_p_read){
      // std::fprintf(stderr, "buffer full, unable to write into buffer, monitor data lose\n");
      return;
    }
    m_vec_ring_ev[next_p_ring_write] = df;
    m_count_ring_write ++;
  }

  std::shared_ptr<altel::TelEvent>& frontBufferEvent(){ //by read thread

    if(m_count_ring_write > m_count_ring_read) {
      uint64_t next_p_ring_read = m_count_ring_read % m_size_ring;
      m_hot_p_read = next_p_ring_read;
      // keep hot read to prevent write-overlapping
      return m_vec_ring_ev[next_p_ring_read];
    }
    else{
      return m_ring_end;
    }
  }

  void popFrontBufferEvent(){ //by read thread
    if(m_count_ring_write > m_count_ring_read) {
      uint64_t next_p_ring_read = m_count_ring_read % m_size_ring;
      m_hot_p_read = next_p_ring_read;
      // keep hot read to prevent write-overlapping
      m_vec_ring_ev[next_p_ring_read].reset();
      m_count_ring_read ++;
    }
  }

  JsonDocument jsd_data_last;

  int drawHook(GLFWwindow* window){
    auto &ev_ref = frontBufferEvent(); //ref only,  no copy, no move
    if(ev_ref){// not nullptr/ring_end
      auto ev = std::move(ev_ref); //moved
      popFrontBufferEvent();

      JsonDocument jsd_data(rapidjson::kObjectType);
      JsonValue js_MeasHits(rapidjson::kArrayType);
      js_MeasHits.Reserve(ev->measHits().size(), jsd_data.GetAllocator());
      for(auto &aMeasHit : ev->measHits()){
        JsonValue js_hit(rapidjson::kArrayType);
        js_hit.Reserve(4, jsd_data.GetAllocator());
        js_hit.PushBack(aMeasHit->u(), jsd_data.GetAllocator());
        js_hit.PushBack(aMeasHit->v(), jsd_data.GetAllocator());
        js_hit.PushBack(aMeasHit->detN(), jsd_data.GetAllocator());
        js_hit.PushBack(1.0, jsd_data.GetAllocator()); // orgin-local-center-mm
        js_MeasHits.PushBack(std::move(js_hit), jsd_data.GetAllocator());
      }
      // JsonUtils::printJsonValue(js_MeasHits, true);
      jsd_data.AddMember("hits", std::move(js_MeasHits), jsd_data.GetAllocator());

      JsonValue js_trajs(rapidjson::kArrayType);
      for(auto &aTraj : ev->TJs){
        JsonValue js_traj(rapidjson::kArrayType);
        js_traj.Reserve(aTraj->trajHits().size(), jsd_data.GetAllocator());
        for(auto &aTrajHit : aTraj->trajHits()){
          JsonValue js_hit(rapidjson::kArrayType);
          js_hit.Reserve(4, jsd_data.GetAllocator());

          js_hit.PushBack(aTrajHit->fitHit()->u(), jsd_data.GetAllocator());
          js_hit.PushBack(aTrajHit->fitHit()->v(), jsd_data.GetAllocator());
          js_hit.PushBack(aTrajHit->detN(), jsd_data.GetAllocator());
          js_hit.PushBack(1.0, jsd_data.GetAllocator()); // 1 orgin-local-center-mm, 0 global

          js_traj.PushBack(std::move(js_hit), jsd_data.GetAllocator());
        }
        js_trajs.PushBack(std::move(js_traj), jsd_data.GetAllocator());
      }
      // JsonUtils::printJsonValue(js_trajs, true);
      jsd_data.AddMember("tracks", std::move(js_trajs), jsd_data.GetAllocator());

      if(jsd_data.HasMember("tracks")){
        telgl->drawTracks(jsd_data);
      }
      if(jsd_data.HasMember("hits")){
        telgl->drawHits(jsd_data);
      }
      if(jsd_data.HasMember("detectors")){
        telgl->drawDetectors(jsd_data);
      }
      else{
        telgl->drawDetectors();
      }
      jsd_data_last = std::move(jsd_data);
    }
    else{
      if(jsd_data_last.HasMember("tracks")){
        telgl->drawTracks(jsd_data_last);
      }
      if(jsd_data_last.HasMember("hits")){
        telgl->drawHits(jsd_data_last);
      }
      if(jsd_data_last.HasMember("detectors")){
        telgl->drawDetectors(jsd_data_last);
      }
      else{
        telgl->drawDetectors();
      }

    }
    return 1;
  }
};
