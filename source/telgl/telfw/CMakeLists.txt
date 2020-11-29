add_library(altel-telfw SHARED src/TelFW.cc)

find_package(glfw3 3.3 REQUIRED)
target_link_libraries(altel-telfw
  PUBLIC glfw altel-telgl
  PRIVATE mycommon
  )

set(LIB_PUBLIC_HEADERS include/TelFW.hh)
set_target_properties(altel-telfw PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

target_include_directories(altel-telfw
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )

install(TARGETS altel-telfw
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

add_subdirectory(exe)
