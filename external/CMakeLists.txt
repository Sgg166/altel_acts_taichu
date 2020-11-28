include(ExternalProject)

if(EXISTS /usr/include/boost169)
  set(BOOST_INCLUDEDIR /usr/include/boost169)
endif()
if(EXISTS /usr/lib64/boost169)
  set(BOOST_LIBRARYDIR /usr/lib64/boost169)
endif()

set(Boost_NO_BOOST_CMAKE ON) # disable new cmake features from Boost 1.70 on
find_package(Boost 1.69 REQUIRED
  COMPONENTS
  system thread filesystem chrono regex program_options unit_test_framework
  # NO_CMAKE_SYSTEM_PATH
  )

# find_package(ROOT REQUIRED HINTS $ENV{ROOTSYS}/cmake NO_CMAKE_SYSTEM_PATH)

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

externalproject_add(
  external_acts
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/acts-core
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}  -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIR} -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIR_RELEASE} -DACTS_BUILD_FATRAS=ON -DACTS_BUILD_ALIGNMENT=ON
  BUILD_ALWAYS true
  BUILD_COMMAND make -j${N_THREADS_BUILD_ACTS}
  INSTALL_COMMAND make install
  COMMAND   echo ""
  COMMAND   echo ""
  COMMAND   echo ">>>>>>>External project ACTS is built and installed<<<<<<<<"
  COMMAND   echo "!!!!Please run cmake again for main project to get updated ACTS installation!!!!"
  COMMAND   echo "eg. cmake .."
  COMMAND   echo ""
  )

set_target_properties(external_acts PROPERTIES EXCLUDE_FROM_ALL TRUE)


# externalproject_add(
#   external_altel_gl
#   SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/altel_gl
#   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
#   BUILD_ALWAYS true
#   BUILD_COMMAND make -j4
#   INSTALL_COMMAND make install
#   )
# set_target_properties(external_altel_gl PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_subdirectory(altel_gl/source/gl)
add_subdirectory(eudaq)
