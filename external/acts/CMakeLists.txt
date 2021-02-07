include(FetchContent)
include(ExternalProject)

cmake_host_system_information(RESULT HOST_MEMORY_MB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT HOST_CPU_CORE QUERY NUMBER_OF_PHYSICAL_CORES)
math(EXPR N_THREADS_BUILD_ACTS "${HOST_MEMORY_MB} / 1000 / 6")
if(N_THREADS_BUILD_ACTS GREATER HOST_CPU_CORE)
  math(EXPR N_THREADS_BUILD_ACTS "${HOST_CPU_CORE} + 1")
endif()

message(STATUS "")
message(STATUS "")
message(STATUS "Host cpu: ${HOST_CPU_CORE} cores")
message(STATUS "Host memory: ${HOST_MEMORY_MB} MB")
message(STATUS "Limited by host memory, ${N_THREADS_BUILD_ACTS} threads can be used to build acts")
message(STATUS "To build ACTS package, use command:  make external_acts")

FetchContent_Declare(
  external_acts
  GIT_REPOSITORY https://github.com/eyiliu/acts-core.git
  GIT_TAG d37a098cf7f0db9933009643ec1921971261cf21 # add_alignment
  # FETCHCONTENT_UPDATES_DISCONNECTED ON
)
FetchContent_Populate(external_acts)

externalproject_add(
  external_acts
  SOURCE_DIR ${external_acts_SOURCE_DIR}
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}  -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIR} -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIR_RELEASE} -DACTS_BUILD_ALIGNMENT=ON
  BUILD_ALWAYS true
  BUILD_COMMAND make -j${N_THREADS_BUILD_ACTS}
  INSTALL_COMMAND make install
  COMMAND   echo ""
  COMMAND   echo ""
  COMMAND   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  COMMAND   echo "!!!External project <ACTS> have been built and installed."
  COMMAND   echo "!!!Please run cmake again for main project to get the updated ACTS installation."
  COMMAND   echo "!!!>> cmake .."
  COMMAND   echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  COMMAND   echo ""
  COMMAND   echo ""
  )

set_target_properties(external_acts PROPERTIES EXCLUDE_FROM_ALL TRUE)
