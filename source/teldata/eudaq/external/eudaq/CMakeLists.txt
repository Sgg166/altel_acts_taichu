
add_definitions("-DEUDAQ_CORE_EXPORTS")
add_definitions("-DEUDAQ_PLATFORM=PF_LINUX")
aux_source_directory(eudaq_src/main/lib/core/src EUDAQ_CORE_SRC)
add_library(eudaq_core SHARED ${EUDAQ_CORE_SRC})

target_link_libraries(eudaq_core
  PUBLIC
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<AND:$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,5.3>>>:stdc++fs>
  # $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
  PRIVATE
  ${CMAKE_DL_LIBS}
  )

target_include_directories(eudaq_core
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/eudaq_src/main/lib/core/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/eudaq_src/main/lib/core/include/eudaq>
  )

# install-tree
install(TARGETS eudaq_core
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

