#pragma once

#include <vector>
#include <memory>

#include <glm/glm.hpp>
#include "gl.h"

#include "myrapidjson.h"

namespace altel{
  class TelGL{
  public:
    struct GeoDataGL{
      GLint   id[4]  {0,0,0,0}; //int, pad
      GLfloat pos[4]  {0,0,0,0}; //vec3, pad
      GLfloat size[4] {30,15,1,0}; //vec3, pad
      GLfloat color[4]{0,1,0,0}; //vec3, pad
      GLfloat pitch[4]{0.028, 0.026, 1, 0}; //pitch x,y, thick z, pad
      GLint   npixel[4]{1024, 512,   1, 0};//pixe number x, y, z/1, pad
      GLfloat miss_alignment[16]{1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1 }; //mat4
    };

    struct TransformDataGL{
      glm::mat4 model;
      glm::mat4 view;
      glm::mat4 proj;
    };

    constexpr static GLuint MAX_ID_SIZE{20};

    ////ubuffer host
    //geometry
    GeoDataGL m_data_geo[MAX_ID_SIZE];
    size_t m_counter_geo{0};
    GLuint m_ubuffer_geo{0};
    GLuint m_bindpoint_geo{0};
    //transform
    TransformDataGL m_data_transform;
    GLuint m_ubuffer_transform{0};
    GLuint m_bindpoint_transform{0};

    ////program telescope
    //shader
    GLuint m_shader_vertex_tel{0};
    GLuint m_shader_geometry_tel{0};
    GLuint m_shader_fragment_tel{0};
    GLuint m_program_tel{0};
    GLuint m_blockindex_geo_tel{0};
    GLuint m_blockindex_transform_tel{0};
    GLint m_location_tel_id{0};
    //vbuffer
    GLuint m_vertex_array_tel{0};
    GLuint m_vbuffer_tel_id{0};
    GLint m_data_tel_id[MAX_ID_SIZE];
    ////

    ////program hit
    //shader
    GLuint m_shader_vertex_hit{0};
    GLuint m_shader_geometry_hit{0};
    GLuint m_shader_fragment_hit{0};
    GLuint m_program_hit{0};
    GLuint m_blockindex_geo_hit{0};
    GLuint m_blockindex_transform_hit{0};
    GLint m_location_hit_pos{0};
    //vbuffer
    GLuint m_vertex_array_hit{0};
    GLuint m_vbuffer_hit_pos{0};
    ////

    ////program track
    //shader
    GLuint m_shader_vertex_track{0};
    GLuint m_shader_fragment_track{0};
    GLuint m_program_track{0};
    GLuint m_blockindex_geo_track{0};
    GLuint m_blockindex_transform_track{0};
    GLint m_location_track_pos{0};
    //vbuffer
    GLuint m_vertex_array_track{0};
    GLuint m_vbuffer_track_pos{0};
    ////

    TelGL(const JsonValue& js);

    ~TelGL();

    void updateGeometry(const JsonValue& js);
    void updateTransform(const JsonValue& js);
    void updateTransform(float cameraX, float cameraY, float cameraZ,
                         float centerX, float centerY, float centerZ,
                         float upvectX, float upvectY, float upvectZ,
                         float povHoriz, float nearDist, float farDist,
                         float ratioWidthHeigh);

    void drawDetectors(const JsonValue& js);
    void drawDetectors();
    void drawHits(const JsonValue& js);
    void drawTracks(const JsonValue& js);
    void draw();


    static JsonDocument createGeometryExample(size_t n=8);
    static JsonDocument createTransformExample();
    static std::string readFile(const std::string& path);
    static JsonDocument createJsonDocument(const std::string& str);
    static JsonValue createJsonValue(JsonAllocator& jsa, const std::string& str);
    static void printJsonValue(const JsonValue& o, bool pretty=false);
  };
}
