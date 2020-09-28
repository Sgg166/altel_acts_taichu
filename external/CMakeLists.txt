include(ExternalProject)

option(ALTEL_BUILD_EXTERNAL_PACKAGE_ACTS "" OFF)
if(ALTEL_BUILD_EXTERNAL_PACKAGE_ACTS)
  message(STATUS "external package acts-core is enabled to build")
  externalproject_add(
    external_acts
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/acts-core
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
    )
endif()
