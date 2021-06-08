add_subdirectory(misc)

aux_source_directory(src LIB_SRC)
add_library(altel-rbcp STATIC ${LIB_SRC} )
target_link_libraries(altel-rbcp
  PUBLIC
  mycommon
  altel-data-event
  PRIVATE
  altel-builtinJSON-headers
  )

target_include_directories(altel-rbcp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  )

install(TARGETS altel-rbcp
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
