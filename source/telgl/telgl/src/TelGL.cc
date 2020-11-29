#include "gl.h"

#include "myrapidjson.h"
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <fstream>
#include <sstream>

#include "TelGL.hh"
#include "glm/ext.hpp"

using namespace altel;

namespace {
  // default tel
  const GLchar* default_vertex_glsl =
#include "TelVertex_glsl.hh"
  ;
  const GLchar* default_geometry_glsl =
#include "TelGeometry_glsl.hh"
  ;
  const GLchar* default_fragment_glsl =
#include "TelFragment_glsl.hh"
  ;

  // default hit
  const GLchar* default_vertex_glsl_hit =
#include "HitVertex_glsl.hh"
  ;
  const GLchar* default_geometry_glsl_hit =
#include "HitGeometry_glsl.hh"
  ;
  const GLchar* default_fragment_glsl_hit =
#include "TelFragment_glsl.hh" // using Tel
  ;

  // default track
  const GLchar* default_vertex_glsl_track =
#include "TrackVertex_glsl.hh"
  ;
  const GLchar* default_fragment_glsl_track =
#include "TelFragment_glsl.hh" // using Tel
  ;

  GLuint createShader(GLenum type, const GLchar* src) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    GLint IsCompiled;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &IsCompiled);
    if(IsCompiled == GL_FALSE){
      GLint maxLength;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);
      std::vector<GLchar> infoLog(maxLength, 0);
      glGetShaderInfoLog(shader, maxLength, &maxLength, infoLog.data());
      std::fprintf(stderr, "ERROR is catched at compile time of opengl GLSL.\n%s", infoLog.data());
      std::fprintf(stderr, "======problematic GLSL code below=====\n%s\n=====end of GLSL code====\n", src);
      throw std::runtime_error("unable to compile glsl");
    }
    return shader;
  }
}

TelGL::TelGL(const JsonValue &js){
  //// create buffer object  (geometry)
  glGenBuffers(1, &m_ubuffer_geo);
  // assign buffer as a vbuffer
  glBindBuffer(GL_UNIFORM_BUFFER, m_ubuffer_geo);
  // assign storage to buffer object
  glNamedBufferData(m_ubuffer_geo, sizeof(m_data_geo), NULL, GL_STATIC_DRAW);
  // assign a bind point to buffer object (and specify the data storage range)
  m_bindpoint_geo = 0;
  glBindBufferRange(GL_UNIFORM_BUFFER, m_bindpoint_geo, m_ubuffer_geo, 0, sizeof(m_data_geo));

  //// create buffer object  (transform)
  glGenBuffers(1, &m_ubuffer_transform);
  // assign buffer as a vbuffer
  glBindBuffer(GL_UNIFORM_BUFFER, m_ubuffer_transform);
  // assign storage to buffer object
  glNamedBufferData(m_ubuffer_transform, sizeof(m_data_transform), NULL, GL_STATIC_DRAW);
  // assign a bind point to buffer object (and specify the data storage range)
  m_bindpoint_transform = 1;
  glBindBufferRange(GL_UNIFORM_BUFFER, m_bindpoint_transform, m_ubuffer_transform, 0, sizeof(m_data_transform));

  //////////////////////////////telescope///////////////////////////////////////
  ////// telescope data
  //// create vertex array object
  glGenVertexArrays(1, &m_vertex_array_tel);
  //// create buffer object
  glGenBuffers(1, &m_vbuffer_tel_id);
  //// assign buffer as a vbuffer in a vertex array object
  glBindVertexArray(m_vertex_array_tel);   // bind vao to GL_ARRAY_BUFFER
  glBindBuffer(GL_ARRAY_BUFFER, m_vbuffer_tel_id);// bind vbo to GL_ARRAY_BUFFER
  //assign storage to buffer object
  glNamedBufferData(m_vbuffer_tel_id, sizeof(m_data_tel_id), NULL, GL_STATIC_DRAW);
  ////// telescope program
  //// shader and variables
  // create program
  m_shader_vertex_tel   = createShader(GL_VERTEX_SHADER,  default_vertex_glsl);
  m_shader_geometry_tel = createShader(GL_GEOMETRY_SHADER, default_geometry_glsl);
  m_shader_fragment_tel = createShader(GL_FRAGMENT_SHADER, default_fragment_glsl);
  m_program_tel = glCreateProgram();
  glAttachShader(m_program_tel, m_shader_vertex_tel);
  glAttachShader(m_program_tel, m_shader_geometry_tel);
  glAttachShader(m_program_tel, m_shader_fragment_tel);
  glLinkProgram(m_program_tel);
  // get attributes from program
  m_blockindex_geo_tel = glGetUniformBlockIndex(m_program_tel, "GeoData");
  m_blockindex_transform_tel = glGetUniformBlockIndex(m_program_tel, "TransformData");
  m_location_tel_id = glGetAttribLocation(m_program_tel, "l");
  glVertexAttribIPointer(m_location_tel_id, 1, GL_INT, sizeof(GLint), 0);
  // connect attributes to buffer
  //// connect uniform block to bind point (ubuffer)
  glUniformBlockBinding(m_program_tel, m_blockindex_geo_tel, m_bindpoint_geo);
  glUniformBlockBinding(m_program_tel, m_blockindex_transform_tel, m_bindpoint_transform);
  //// connect vertex attribute to vertex array object (vbuffer)
  glBindVertexArray(m_vertex_array_tel);
  glEnableVertexAttribArray(m_location_tel_id); //use currently bound vertex array object for the operation
  //////////////////////////////end of telescope////////////////////////////////////

  //////////////////////////////hit///////////////////////////////////////
  ////// hit data
  //// create vertex array object
  glGenVertexArrays(1, &m_vertex_array_hit);
  //// create buffer object
  glGenBuffers(1, &m_vbuffer_hit_pos);
  //// assign buffer as a vbuffer in a vertex array object
  glBindVertexArray(m_vertex_array_hit);   // bind vao to GL_ARRAY_BUFFER
  glBindBuffer(GL_ARRAY_BUFFER, m_vbuffer_hit_pos);// bind vbo to GL_ARRAY_BUFFER
  //assign storage to buffer object, MAX 1024 hits
  glNamedBufferData(m_vbuffer_hit_pos, sizeof(GLfloat)*4*1024, NULL, GL_STATIC_DRAW);
  ////// program
  //// shader and variables
  // create program
  m_shader_vertex_hit   = createShader(GL_VERTEX_SHADER,   default_vertex_glsl_hit);
  m_shader_geometry_hit = createShader(GL_GEOMETRY_SHADER, default_geometry_glsl_hit);
  m_shader_fragment_hit = createShader(GL_FRAGMENT_SHADER, default_fragment_glsl_hit);
  m_program_hit = glCreateProgram();
  glAttachShader(m_program_hit, m_shader_vertex_hit);
  glAttachShader(m_program_hit, m_shader_geometry_hit);
  glAttachShader(m_program_hit, m_shader_fragment_hit);
  glLinkProgram(m_program_hit);
  // get attributes from program
  m_blockindex_geo_hit = glGetUniformBlockIndex(m_program_hit, "GeoData");
  m_blockindex_transform_hit = glGetUniformBlockIndex(m_program_hit, "TransformData");
  m_location_hit_pos = glGetAttribLocation(m_program_hit, "pos");
  glVertexAttribPointer(m_location_hit_pos, 4, GL_FLOAT, GL_FALSE, sizeof(GLint)*4, 0);
  // connect attributes to buffer
  //// connect uniform block to bind point (ubuffer)
  glUniformBlockBinding(m_program_hit, m_blockindex_geo_hit, m_bindpoint_geo);
  glUniformBlockBinding(m_program_hit, m_blockindex_transform_hit, m_bindpoint_transform);
  //// connect vertex attribute to vertex array object (vbuffer)
  glBindVertexArray(m_vertex_array_hit);
  glEnableVertexAttribArray(m_location_hit_pos); //use currently bound vertex array object for the operation
  //////////////////////////////end of hit////////////////////////////////////


  //////////////////////////////track///////////////////////////////////////
  ////// data
  //// create vertex array object
  glGenVertexArrays(1, &m_vertex_array_track);
  //// create buffer object
  glGenBuffers(1, &m_vbuffer_track_pos);
  //// assign buffer as a vbuffer in a vertex array object
  glBindVertexArray(m_vertex_array_track);   // bind vao to GL_ARRAY_BUFFER
  glBindBuffer(GL_ARRAY_BUFFER, m_vbuffer_track_pos);// bind vbo to GL_ARRAY_BUFFER
  //assign storage to buffer object
  glNamedBufferData(m_vbuffer_track_pos, sizeof(GLfloat)*4*1024, NULL, GL_STATIC_DRAW); //TODO
  ////// program
  //// shader and variables
  // create program
  m_shader_vertex_track   = createShader(GL_VERTEX_SHADER,   default_vertex_glsl_track);
  m_shader_fragment_track = createShader(GL_FRAGMENT_SHADER, default_fragment_glsl_track);
  m_program_track = glCreateProgram();
  glAttachShader(m_program_track, m_shader_vertex_track);
  glAttachShader(m_program_track, m_shader_fragment_track);
  glLinkProgram(m_program_track);
  // get attributes from program
  m_blockindex_geo_track = glGetUniformBlockIndex(m_program_track, "GeoData");
  m_blockindex_transform_track = glGetUniformBlockIndex(m_program_track, "TransformData");
  m_location_track_pos = glGetAttribLocation(m_program_track, "pos");
  glVertexAttribPointer(m_location_track_pos, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*4, 0);  //TODO
  // connect attributes to buffer
  //// connect uniform block to bind point (ubuffer)
  glUniformBlockBinding(m_program_track, m_blockindex_geo_track, m_bindpoint_geo);
  glUniformBlockBinding(m_program_track, m_blockindex_transform_track, m_bindpoint_transform);
  //// connect vertex attribute to vertex array object (vbuffer)
  glBindVertexArray(m_vertex_array_track);
  glEnableVertexAttribArray(m_location_track_pos); //use currently bound vertex array object for the operation
  //////////////////////////////end of track////////////////////////////////////

  if(!js.IsNull()){
    updateGeometry(js);
  }

}

TelGL::~TelGL(){
  ////// telescope
  if(m_program_tel) {glDeleteProgram(m_program_tel); m_program_tel = 0;}
  if(m_shader_fragment_tel){glDeleteShader(m_shader_fragment_tel); m_shader_fragment_tel = 0;}
  if(m_shader_geometry_tel){glDeleteShader(m_shader_geometry_tel); m_shader_geometry_tel = 0;}
  if(m_shader_vertex_tel){glDeleteShader(m_shader_vertex_tel); m_shader_vertex_tel = 0;}
  if(m_vertex_array_tel){glDeleteVertexArrays(1, &m_vertex_array_tel); m_vertex_array_tel = 0;}
  if(m_vbuffer_tel_id){glDeleteBuffers(1, &m_vbuffer_tel_id); m_vbuffer_tel_id = 0;}
  //
  ////// hit
  if(m_program_hit) {glDeleteProgram(m_program_hit); m_program_hit = 0;}
  if(m_shader_fragment_hit){glDeleteShader(m_shader_fragment_hit); m_shader_fragment_hit = 0;}
  if(m_shader_geometry_hit){glDeleteShader(m_shader_geometry_hit); m_shader_geometry_hit = 0;}
  if(m_shader_vertex_hit){glDeleteShader(m_shader_vertex_hit); m_shader_vertex_hit = 0;}
  if(m_vertex_array_hit){glDeleteVertexArrays(1, &m_vertex_array_hit); m_vertex_array_hit = 0;}
  if(m_vbuffer_hit_pos){glDeleteBuffers(1, &m_vbuffer_hit_pos); m_vbuffer_hit_pos = 0;}
  //

  ////// track
  if(m_program_track) {glDeleteProgram(m_program_track); m_program_track = 0;}
  if(m_shader_fragment_track){glDeleteShader(m_shader_fragment_track); m_shader_fragment_track = 0;}
  if(m_shader_vertex_track){glDeleteShader(m_shader_vertex_track); m_shader_vertex_track = 0;}
  if(m_vertex_array_track){glDeleteVertexArrays(1, &m_vertex_array_track); m_vertex_array_track = 0;}
  if(m_vbuffer_track_pos){glDeleteBuffers(1, &m_vbuffer_track_pos); m_vbuffer_track_pos = 0;}
  //

  // ubuffer
}

void TelGL::updateTransform(float cameraX, float cameraY, float cameraZ,
                            float centerX, float centerY, float centerZ,
                            float upvectX, float upvectY, float upvectZ,
                            float fovHoriz, float nearDist, float farDist,
                            float ratioWidthHigh){
  m_data_transform.model = glm::mat4( 1.0f );
  m_data_transform.view = glm::lookAt(glm::vec3(cameraX, cameraY, cameraZ),    // eye position
                                      glm::vec3(centerX, centerY, centerZ),    // object center
                                      glm::vec3(upvectX, upvectY, upvectZ));   // up vector
  m_data_transform.proj = glm::perspective(glm::radians(fovHoriz), ratioWidthHigh, nearDist, farDist);

  glNamedBufferSubData(m_ubuffer_transform, 0, sizeof(m_data_transform), &m_data_transform);
}


void TelGL::updateTransform(const JsonValue& js){
  if(!js.HasMember("trans")){
    std::fprintf(stderr, "unable to find \"trans\" key from JS\n");
    printJsonValue(js, true);
    throw;
    return;
  }
  const auto &js_trans = js["trans"];
  glm::mat4 mvp[3]; // 0 model, 1 view, 2 proj

  // const auto &js_model = js_trans["model"];
  mvp[0] = glm::mat4( 1.0f );

  const auto &js_lookat = js_trans["lookat"];
  const auto &js_eye = js_lookat["eye"];
  const auto &js_center = js_lookat["center"];
  const auto &js_up = js_lookat["up"];
  mvp[1] = glm::lookAt(glm::vec3(js_eye["x"].GetDouble(), js_eye["y"].GetDouble(), js_eye["z"].GetDouble()),    // eye position
                       glm::vec3(js_center["x"].GetDouble(), js_center["y"].GetDouble(), js_center["z"].GetDouble()), // object center
                       glm::vec3(js_up["x"].GetDouble(), js_up["y"].GetDouble(), js_up["z"].GetDouble()));   // up vector

  const auto &js_pers = js_trans["persp"];
  const auto &js_fov = js_pers["fov"];
  const auto &js_ratio = js_pers["ratio"];
  const auto &js_near = js_pers["near"];
  const auto &js_far = js_pers["far"];
  mvp[2] = glm::perspective(glm::radians(js_fov.GetDouble()), js_ratio.GetDouble(), js_near.GetDouble(), js_far.GetDouble());

  glNamedBufferSubData(m_ubuffer_transform, 0, sizeof(mvp), mvp);
}


void TelGL::updateGeometry(const JsonValue& js){
  if(!js.HasMember("geometry")){
    std::fprintf(stderr, "unable to find \"geomerty\" key for detector geomerty from JS\n");
    return;
  }

  const auto &js_geo = js["geometry"];
  const auto &js_dets = js_geo["detectors"];
  size_t n = 0;
  for(const auto& js_det: js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cx = js_det["center"]["x"].GetDouble();
    double cy = js_det["center"]["y"].GetDouble();
    double cz = js_det["center"]["z"].GetDouble();
    double rx = js_det["rotation"]["x"].GetDouble();
    double ry = js_det["rotation"]["y"].GetDouble();
    double rz = js_det["rotation"]["z"].GetDouble();
    double ptx = js_det["pitch"]["x"].GetDouble();
    double pty = js_det["pitch"]["y"].GetDouble();
    double ptz = js_det["pitch"]["z"].GetDouble();
    double px = js_det["pixel"]["x"].GetDouble();
    double py = js_det["pixel"]["y"].GetDouble();
    double pz = js_det["pixel"]["z"].GetDouble();
    double sx = js_det["size"]["x"].GetDouble();
    double sy = js_det["size"]["y"].GetDouble();
    double sz = js_det["size"]["z"].GetDouble();

    double colorR = 0;
    double colorG = 0;
    double colorB = 1.0;

    if(js_det.HasMember("color")){
      double colorR = js_det["color"]["r"].GetDouble();
      double colorG = js_det["color"]["g"].GetDouble();
      double colorB = js_det["color"]["b"].GetDouble();
    }

    m_data_tel_id[n] = id; //
    m_data_geo[n].id[0] = id;
    m_data_geo[n].pos[0]=  cx;
    m_data_geo[n].pos[1]=  cy;
    m_data_geo[n].pos[2]=  cz;
    m_data_geo[n].pos[3]=  0;
    m_data_geo[n].size[0] = sx;
    m_data_geo[n].size[1] = sy;
    m_data_geo[n].size[2] = sz;
    m_data_geo[n].size[3] = 0;
    m_data_geo[n].color[0]= colorR;
    m_data_geo[n].color[1]= colorG;
    m_data_geo[n].color[2]= colorB;
    m_data_geo[n].color[3]= 0;
    m_data_geo[n].pitch[0]= ptx;
    m_data_geo[n].pitch[1]= pty;
    m_data_geo[n].pitch[2]= ptz;
    m_data_geo[n].pitch[3]= 0;
    m_data_geo[n].npixel[0]= px;
    m_data_geo[n].npixel[1]= py;
    m_data_geo[n].npixel[2]= pz;
    m_data_geo[n].npixel[3]= 0;
    n++;
    if(n>=sizeof(m_data_geo)/sizeof(TelGL::GeoDataGL)){
      std::fprintf(stdout, "TelGL:: Max size of geo is reached, ignoring later detectors!\n");
      break;
    }
  }
  m_counter_geo = n;
  glNamedBufferSubData(m_ubuffer_geo, 0, sizeof(TelGL::GeoDataGL)*m_counter_geo, &m_data_geo[0]);
  //std::fprintf(stdout, "set up geometry data with %i layers\n", m_counter_geo);
}

void TelGL::drawDetectors(const JsonValue& js){
  if(js.HasMember("detectors")){
    const auto &js_data = js["detectors"];
    if(!js_data.IsArray()){
      std::fprintf(stderr, "unable to to decode detectors from JS\n");
      return;
    }
    glUseProgram(m_program_tel);
    glBindVertexArray(m_vertex_array_tel);
    std::vector<GLint> gldata(js_data.Size());
    auto it = gldata.begin();
    for(const auto &js_det : js_data.GetArray()){
      *(it++) = js_det.GetInt();
    }
    glNamedBufferSubData(m_vbuffer_tel_id, 0, sizeof(GLint)*gldata.size(), gldata.data());
    glDrawArrays(GL_POINTS, 0, gldata.size());
  }
}

void TelGL::drawDetectors(){
  glUseProgram(m_program_tel);
  glBindVertexArray(m_vertex_array_tel);
  glNamedBufferSubData(m_vbuffer_tel_id, 0, sizeof(GLint)*m_counter_geo, m_data_tel_id);
  glDrawArrays(GL_POINTS, 0, m_counter_geo);
}

// TODO: change glsl from ivec3 to vec4
void TelGL::drawHits(const JsonValue& js){
  if(js.HasMember("hits")){
    const auto &js_data = js["hits"];
    if(!js_data.IsArray()){
      std::fprintf(stderr, "unable to to decode hits from JS\n");
      printJsonValue(js, true);
      return;
    }
    glUseProgram(m_program_hit);
    glBindVertexArray(m_vertex_array_hit);
    // TODO: glsl is vec3
    std::vector<GLfloat> gldata(4*js_data.Size());
    auto it = gldata.begin();
    for(const auto &js_hit : js_data.GetArray()){
      *(it++) = js_hit[0].GetDouble();
      *(it++) = js_hit[1].GetDouble();
      *(it++) = js_hit[2].GetDouble();
      *(it++) = js_hit[3].GetDouble();
    }
    glNamedBufferSubData(m_vbuffer_hit_pos, 0, sizeof(GLfloat)*gldata.size(), gldata.data());
    glDrawArrays(GL_POINTS, 0, gldata.size()/4);
  }
}

void TelGL::drawTracks(const JsonValue& js){
  if(js.HasMember("tracks")){
    const auto &js_data = js["tracks"];
    if(!js_data.IsArray()){
      std::fprintf(stderr, "unable to to decode tracks from JS\n");
      printJsonValue(js, true);
      return;
    }
    glUseProgram(m_program_track);
    glBindVertexArray(m_vertex_array_track);
    for(const auto &js_track : js_data.GetArray()){
      std::vector<GLfloat> gldata(4*js_track.Size());
      auto it = gldata.begin();
      for(const auto &js_hit : js_track.GetArray()){
        *(it++) = js_hit[0].GetDouble();
        *(it++) = js_hit[1].GetDouble();
        *(it++) = js_hit[2].GetDouble();
        *(it++) = js_hit[3].GetDouble();
      }
      glNamedBufferSubData(m_vbuffer_track_pos, 0, sizeof(GLfloat)*gldata.size(), gldata.data());
      glDrawArrays(GL_LINE_STRIP, 0, gldata.size()/4);
    }
  }
}


JsonDocument TelGL::createGeometryExample(size_t num){
  static constexpr const float colors[][3]=
    {
     {1.0, 0.0, 0.0}, //red
     {0.0, 1.0, 0.0}, //green
     {0.0, 0.0, 1.0}, //blue
     {0.0, 0.5, 0.5}, //darkcyan
     {1.0, 0.0, 1.0}, //magenta
     {0.5, 0.5, 0.0}, //darkyellow
     {0.5, 0.0, 0.0}, //darkred
     {0.0, 0.5, 0.0}, //darkgreen
     {0.0, 0.0, 0.5}, //darkblue
     {0.5, 0.0, 0.5}, //darkmagenta
     {0.0, 0.0, 0.0}, //black
     {0.0, 1.0, 1.0}, //cyan
     {1.0, 1.0, 0.0}}; //yellow
  static constexpr const size_t ncolor=
    sizeof(colors)/sizeof(float)/3;

  JsonDocument jsdoc;
  jsdoc.SetObject();
  JsonAllocator& jsa= jsdoc.GetAllocator();
  using rapidjson::kObjectType;
  using rapidjson::kArrayType;
  JsonValue js_geo(kObjectType);
  JsonValue js_dets(kArrayType);
  for(size_t i=0; i<num; i++){
    JsonValue js_det(kObjectType);
    js_det.AddMember("id", i, jsa);

    JsonValue js_size(kObjectType);
    js_size.AddMember("x", 0.02924*1024., jsa);
    js_size.AddMember("y", 0.02688*512., jsa);
    js_size.AddMember("z", 1., jsa);
    js_det.AddMember("size", std::move(js_size), jsa);

    JsonValue js_pitch(kObjectType);
    js_pitch.AddMember("x", 0.02924, jsa);
    js_pitch.AddMember("y", 0.02688, jsa);
    js_pitch.AddMember("z", 1., jsa);
    js_det.AddMember("pitch", std::move(js_pitch), jsa);

    JsonValue js_pixel(kObjectType);
    js_pixel.AddMember("x", 1024, jsa);
    js_pixel.AddMember("y", 512, jsa);
    js_pixel.AddMember("z", 1, jsa);
    js_det.AddMember("pixel", std::move(js_pixel), jsa);

    JsonValue js_color(kObjectType);
    js_color.AddMember("r", colors[i%ncolor][0], jsa);
    js_color.AddMember("g", colors[i%ncolor][1], jsa);
    js_color.AddMember("b", colors[i%ncolor][2], jsa);
    js_det.AddMember("color", std::move(js_color), jsa);

    JsonValue js_center(kObjectType);
    js_center.AddMember("x", 0., jsa);
    js_center.AddMember("y", 0., jsa);
    js_center.AddMember("z", 40.*i, jsa);
    js_det.AddMember("center", std::move(js_center), jsa);

    JsonValue js_rotation(kObjectType);
    js_rotation.AddMember("x", 0., jsa);
    js_rotation.AddMember("y", 0., jsa);
    js_rotation.AddMember("z", 0., jsa);
    js_det.AddMember("rotation", std::move(js_rotation), jsa);

    js_dets.PushBack(js_det, jsa);
  }
  js_geo.AddMember("detectors", std::move(js_dets), jsa);
  jsdoc.AddMember("geometry", std::move(js_geo), jsa);

  return jsdoc;
}


JsonDocument TelGL::createTransformExample(){
  JsonDocument jsdoc_trans;
  jsdoc_trans.SetObject();
  JsonAllocator& jsa= jsdoc_trans.GetAllocator();

  using rapidjson::kObjectType;
  JsonValue js_trans(kObjectType);

  JsonValue js_lookat(kObjectType);
  JsonValue js_eye(kObjectType);
  js_eye.AddMember("x", 0., jsa);
  js_eye.AddMember("y", 0., jsa);
  js_eye.AddMember("z", -1000., jsa);
  js_lookat.AddMember("eye", std::move(js_eye), jsa);
  JsonValue js_center(kObjectType);
  js_center.AddMember("x", 0., jsa);
  js_center.AddMember("y", 0., jsa);
  js_center.AddMember("z", 0., jsa);
  js_lookat.AddMember("center", std::move(js_center), jsa);
  JsonValue js_up(kObjectType);
  js_up.AddMember("x", 0., jsa);
  js_up.AddMember("y", 1., jsa);
  js_up.AddMember("z", 0., jsa);
  js_lookat.AddMember("up", std::move(js_up), jsa);
  js_trans.AddMember("lookat", std::move(js_lookat), jsa);

  JsonValue js_persp(kObjectType);
  js_persp.AddMember("fov", 60., jsa);
  js_persp.AddMember("ratio",2., jsa);
  js_persp.AddMember("near",0.1, jsa);
  js_persp.AddMember("far", 2000., jsa);
  js_trans.AddMember("persp", std::move(js_persp), jsa);

  jsdoc_trans.AddMember("trans", std::move(js_trans), jsa);
  return jsdoc_trans;
}

JsonDocument TelGL::createJsonDocument(const std::string& str){
  std::istringstream iss(str);
  rapidjson::IStreamWrapper isw(iss);
  JsonDocument jsdoc;
  jsdoc.ParseStream(isw);
  if(jsdoc.HasParseError()){
    std::fprintf(stderr, "rapidjson error<%s> when parsing input data at offset %llu\n",
                 rapidjson::GetParseError_En(jsdoc.GetParseError()), jsdoc.GetErrorOffset());
    throw std::runtime_error("TelGL::createJsonDocument has ParseError\n");
  }
  return jsdoc;
}

JsonValue TelGL::createJsonValue(JsonAllocator& jsa, const std::string& str){
  JsonDocument jsdoc = createJsonDocument(str);
  JsonValue js;
  js.CopyFrom(jsdoc, jsa);
  return js;
}

std::string TelGL::readFile(const std::string& path){
  std::ifstream ifs(path);
  if(!ifs.good()){
    std::fprintf(stderr, "unable to raed file path<%s>\n", path);
    throw std::runtime_error("File open error\n");
  }
  std::string str;
  str.assign((std::istreambuf_iterator<char>(ifs) ),
             (std::istreambuf_iterator<char>()));
  return str;
}

void TelGL::printJsonValue(const JsonValue& o, bool pretty){
  rapidjson::StringBuffer sb;
  if(pretty){
    rapidjson::PrettyWriter<rapidjson::StringBuffer> w(sb);
    o.Accept(w);
  }
  else{
    rapidjson::Writer<rapidjson::StringBuffer> w(sb);
    o.Accept(w);
  }
  rapidjson::PutN(sb, '\n', 1);
  std::fwrite(sb.GetString(), 1, sb.GetSize(), stdout);
}


