add_subdirectory(event)
add_subdirectory(root)
# add_subdirectory(json)

add_library(altel-data INTERFACE)

target_link_libraries(altel-data INTERFACE
  mycommon
  altel-data-event
  altel-data-root
  )

install(TARGETS altel-data
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
