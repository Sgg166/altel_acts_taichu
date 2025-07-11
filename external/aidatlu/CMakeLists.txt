include(FetchContent)
include(ExternalProject)
FetchContent_Declare(
  external_aidatlu
  GIT_REPOSITORY https://github.com/eyiliu/aidatlu.git
  GIT_TAG 43d9376a3642555ffb972067e61d5be28593a709
  # FETCHCONTENT_UPDATES_DISCONNECTED ON
)
FetchContent_Populate(external_aidatlu)

externalproject_add(
  external_aidatlu
  SOURCE_DIR ${external_aidatlu_SOURCE_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  BUILD_ALWAYS true
  BUILD_COMMAND make -j2
  INSTALL_COMMAND make install
  COMMAND   echo ""
  COMMAND   echo ""
  )

set_target_properties(external_aidatlu PROPERTIES EXCLUDE_FROM_ALL TRUE)
