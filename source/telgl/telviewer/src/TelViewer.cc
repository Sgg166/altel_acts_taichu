#include "TelGL.hh"
#include "TelViewer.hh"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

using namespace altel;

TelViewer::~TelViewer(){
  stopAsyncLoop();
}

void TelViewer::clearObjects(){ // by write thread
  uint64_t count_ring_read = m_count_ring_read;
  m_count_ring_write = count_ring_read;
}

void TelViewer::pushObject(std::shared_ptr<JsonValue> dp){ //by write thread
  uint64_t next_p_ring_write = m_count_ring_write % m_size_ring;
  if(next_p_ring_write == m_hot_p_read){
    std::fprintf(stderr, "buffer full, unable to write into buffer, monitor data lose\n");
    return;
  }
  m_vec_ring_ev[next_p_ring_write] = dp;
  m_count_ring_write ++;
}

std::shared_ptr<JsonValue>& TelViewer::frontObject(){ //by read thread
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

void TelViewer::popFrontObject(){ //by read thread
  if(m_count_ring_write > m_count_ring_read) {
    uint64_t next_p_ring_read = m_count_ring_read % m_size_ring;
    m_hot_p_read = next_p_ring_read;
    // keep hot read to prevent write-overlapping
    m_vec_ring_ev[next_p_ring_read].reset();
    m_count_ring_read ++;
  }
}

void TelViewer::stopAsyncLoop(){
  m_is_async_thread_running = false;
  if(m_fut_async_thread.valid())
    m_fut_async_thread.get();
}

void TelViewer::startAsyncLoop(){
  if(m_is_async_thread_running){
    std::fprintf(stderr, "unable to start up new async thread.  already running\n");
  }
  m_vec_ring_ev.clear();
  m_vec_ring_ev.resize(m_size_ring);
  m_count_ring_write = 0;
  m_count_ring_read = 0;
  m_hot_p_read = m_size_ring -1; // tail
  m_fut_async_thread = std::async(std::launch::async, &TelViewer::asyncLoop, this);
}

uint64_t TelViewer::asyncLoop(){
  glfwSetErrorCallback([](int error, const char* description){
                         fprintf(stderr, "Error: %s\n", description);});
  if (!glfwInit()){
    std::fprintf(stderr, "unable to init by glfw\n");
    return 0;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_SAMPLES, 4);

  GLFWwindow* window = glfwCreateWindow((int)800, (int)600, "OpenGL telescope example", NULL, NULL);
  if (!window){
    glfwTerminate();
    std::fprintf(stderr, "unable to create window by glfw\n");
    return 0;
  }
  // now no error return, set flag to true;
  bool &flag_running = m_is_async_thread_running;
  flag_running = true;

  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height){
                                           glViewport(0, 0, width, height);});


  //TODO: create geo js
  JsonValue js_geo;
  TelGL telgl(js_geo);
  uint64_t n_gl_frame = 0;

  float deltaTime = 0.0f; // time between current frame and last frame
  float lastFrame = 0.0f;
  glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  400.0f);
  glm::vec3 worldCenter = glm::vec3(0.0f, 0.0f,  0.0f);
  glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);
  JsonDocument jsconf_trans = altel::TelGL::createTransformExample();
  while(flag_running && !glfwWindowShouldClose(window)){
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
      glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
      cameraPos -= (cameraPos - worldCenter) * (deltaTime * 0.1f);
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
      cameraPos += (cameraPos - worldCenter) * (deltaTime * 0.1f);
    //note: rotation does not work as expectation, BUT WORKS SOMEHOW
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
      cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
      cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(0.0f, 1.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
      cameraPos -= glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
      cameraPos += glm::cross((cameraPos - worldCenter), glm::vec3(1.0f, 0.0f,  0.0f)) * (deltaTime * 2.0f);

    // pausing
    if(glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS){
      continue;
    }

    // continue; // debugging
    auto &js_ref = frontObject(); //ref only,  no copy, no move
    if(js_ref == m_ring_end){// not nullptr/ring_end
      continue;
    }
    auto jssp = std::move(js_ref); //moved
    popFrontObject();

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    const float ratio = width / (float) height;
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    jsconf_trans["trans"]["lookat"]["eye"]["x"]= cameraPos[0];
    jsconf_trans["trans"]["lookat"]["eye"]["y"]= cameraPos[1];
    jsconf_trans["trans"]["lookat"]["eye"]["z"]= cameraPos[2];
    jsconf_trans["trans"]["lookat"]["center"]["x"]= worldCenter[0];
    jsconf_trans["trans"]["lookat"]["center"]["y"]= worldCenter[1];
    jsconf_trans["trans"]["lookat"]["center"]["z"]= worldCenter[2];
    jsconf_trans["trans"]["lookat"]["up"]["x"]= cameraUp[0];
    jsconf_trans["trans"]["lookat"]["up"]["y"]= cameraUp[1];
    jsconf_trans["trans"]["lookat"]["up"]["z"]= cameraUp[2];
    jsconf_trans["trans"]["persp"]["fov"]= 60;
    jsconf_trans["trans"]["persp"]["ratio"]= ratio;
    jsconf_trans["trans"]["persp"]["near"]= 0.1;
    jsconf_trans["trans"]["persp"]["far"]= 2000;

    telgl.updateTransform(jsconf_trans);

    const JsonValue& js_data = *jssp;
    if(js_data.HasMember("tracks")){
      telgl.drawTracks(js_data);
    }
    if(js_data.HasMember("hits")){
      telgl.drawHits(js_data);
    }
    if(js_data.HasMember("detectors")){
      telgl.drawDetectors(js_data);
    }
    else{
      telgl.drawDetectors();
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
    n_gl_frame ++;
  }

  glfwDestroyWindow(window);
  glfwTerminate();
  flag_running = false;
  return n_gl_frame;
}
