unset(CURRENT_TARGETS_LIST)

add_executable(altelMilleBin_json TelMille_main.cpp)
target_link_libraries(altelMilleBin_json PRIVATE mycommon altel-mille)
list(APPEND CURRENT_TARGETS_LIST altelMilleBin_json)

add_executable(altelGeoUpdate GeometryFileUpdate_main.cpp)
target_link_libraries(altelGeoUpdate PRIVATE mycommon altel-mille)
list(APPEND CURRENT_TARGETS_LIST altelGeoUpdate)

install(TARGETS ${CURRENT_TARGETS_LIST}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
