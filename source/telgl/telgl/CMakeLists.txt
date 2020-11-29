add_library(altel-telgl SHARED src/TelGL.cc)

set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)

target_link_libraries(altel-telgl
  PUBLIC galogen  OpenGL::GL
  PRIVATE altel-glsl mycommon
  )

set(LIB_PUBLIC_HEADERS include/TelGL.hh)
set_target_properties(altel-telgl PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

target_include_directories(altel-telgl
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )

install(TARGETS altel-telgl
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
