add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/altel_tele_conf_json.hh
  COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/cpprawstring.cmake ${CMAKE_CURRENT_BINARY_DIR}/altel_tele_conf_json.hh ${CMAKE_CURRENT_SOURCE_DIR}/altel_tele_conf.json
  )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/altel_layer_conf_json.hh
  COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/cpprawstring.cmake ${CMAKE_CURRENT_BINARY_DIR}/altel_layer_conf_json.hh ${CMAKE_CURRENT_SOURCE_DIR}/altel_layer_conf.json
  )

add_custom_target(altel-builtinJSON-custom ALL DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/altel_tele_conf_json.hh
  ${CMAKE_CURRENT_BINARY_DIR}/altel_layer_conf_json.hh)

add_library(altel-builtinJSON-headers INTERFACE)

add_dependencies(altel-builtinJSON-headers altel-builtinJSON-custom)
target_include_directories(altel-builtinJSON-headers INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

install(TARGETS altel-builtinJSON-headers
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
