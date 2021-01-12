include(FetchContent)

FetchContent_Declare(
  external_eudaq_core_only
  GIT_REPOSITORY https://github.com/eudaq/eudaq.git
  GIT_TAG fc77bf805509e910091917bfa4072ced3f5b3480 # 2.4.4 + 1 commit
  # FETCHCONTENT_UPDATES_DISCONNECTED ON
)

FetchContent_Populate(external_eudaq_core_only)

add_definitions("-DEUDAQ_CORE_EXPORTS")
add_definitions("-DEUDAQ_PLATFORM=PF_LINUX")

aux_source_directory(${external_eudaq_core_only_SOURCE_DIR}/main/lib/core/src EUDAQ_CORE_SRC)
add_library(eudaq_core_only SHARED ${EUDAQ_CORE_SRC})

target_link_libraries(eudaq_core_only
  PUBLIC
  $<$<AND:$<CXX_COMPILER_ID:GNU>,$<AND:$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>,$<VERSION_GREATER:$<CXX_COMPILER_VERSION>,5.3>>>:stdc++fs>
  PRIVATE
  ${CMAKE_DL_LIBS}
  )

target_include_directories(eudaq_core_only
  PUBLIC
  $<BUILD_INTERFACE:${external_eudaq_core_only_SOURCE_DIR}/main/lib/core/include>
  $<BUILD_INTERFACE:${external_eudaq_core_only_SOURCE_DIR}/main/lib/core/include/eudaq>
  )

# install-tree
install(TARGETS eudaq_core_only
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

