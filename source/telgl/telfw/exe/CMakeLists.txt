find_package(glfw3 3.3 REQUIRED)

add_executable(TelGL_test TelGL_test.cxx)
target_link_libraries(TelGL_test PRIVATE mycommon altel-telgl altel-telfw glfw)
set_target_properties(TelGL_test PROPERTIES CXX_STANDARD 17)


add_executable(TelFW_test TelFW_test.cxx)
target_link_libraries(TelFW_test PRIVATE mycommon altel-telgl altel-telfw glfw)
set_target_properties(TelFW_test PROPERTIES CXX_STANDARD 17)


install(TARGETS TelGL_test TelFW_test
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
