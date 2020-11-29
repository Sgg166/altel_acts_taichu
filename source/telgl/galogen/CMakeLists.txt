set(LIB_SRC gl.c)
add_library(galogen STATIC ${LIB_SRC})

# target_link_libraries(galogen PRIVATE GL) 
target_include_directories(galogen
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  )

set_source_files_properties(gl.c PROPERTIES LANGUAGE CXX)

set(LIB_PUBLIC_HEADERS gl.h)
set_target_properties(galogen PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")
set_target_properties(galogen PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS galogen
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
