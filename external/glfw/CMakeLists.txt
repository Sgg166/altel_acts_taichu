include(FetchContent)
include(ExternalProject)
FetchContent_Declare(
  external_glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG 3.3.2
  # FETCHCONTENT_UPDATES_DISCONNECTED ON
)
FetchContent_Populate(external_glfw)

externalproject_add(
  external_glfw
  SOURCE_DIR ${external_glfw_SOURCE_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  BUILD_ALWAYS true
  BUILD_COMMAND make -j2
  INSTALL_COMMAND make install
  COMMAND   echo ""
  COMMAND   echo ""
  )

set_target_properties(external_glfw PROPERTIES EXCLUDE_FROM_ALL TRUE)
