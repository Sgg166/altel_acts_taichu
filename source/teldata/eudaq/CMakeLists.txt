add_subdirectory(external)

aux_source_directory(src LIB_SRC)
add_library(altel-data-eudaq SHARED ${LIB_SRC})

target_include_directories( altel-data-eudaq
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )

target_link_libraries(
  altel-data-eudaq 
  PUBLIC altel-data-event eudaq_core_only
  PRIVATE mycommon
  )

# install-tree
install(TARGETS altel-data-eudaq
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
