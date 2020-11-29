add_library(altel-telviewer SHARED src/TelViewer.cc)

find_package(glfw3 3.3 REQUIRED)
target_link_libraries(altel-telviewer
  PUBLIC 
  PRIVATE mycommon glfw altel-telgl
  )

set(LIB_PUBLIC_HEADERS include/TelViewer.hh)
set_target_properties(altel-telviewer PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

target_include_directories(altel-telviewer
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )

install(TARGETS altel-telviewer
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
