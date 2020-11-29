aux_source_directory(src LIB_SRC)
add_library(altel-acts SHARED ${LIB_SRC})

target_include_directories( altel-acts
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  # $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )

target_link_libraries(
  altel-acts 
  PUBLIC ActsCore altel-data
  PRIVATE mycommon
  )

file(GLOB LIB_ALL_HEADERS "include/*.hpp")
set(LIB_PUBLIC_HEADERS "${LIB_ALL_HEADERS}" )
set_target_properties(altel-acts PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")


# install-tree
install(TARGETS altel-acts
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

