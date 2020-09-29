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

find_package(ROOT REQUIRED HINTS $ENV{ROOTSYS}/cmake NO_CMAKE_SYSTEM_PATH)

externalproject_add(
  external_acts
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/acts-core
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}  -DBOOST_INCLUDEDIR=${Boost_INCLUDE_DIR} -DBOOST_LIBRARYDIR=${Boost_LIBRARY_DIR_RELEASE} -DACTS_BUILD_FATRAS=ON -DACTS_BUILD_DIGITIZATION_PLUGIN=ON -DACTS_BUILD_IDENTIFICATION_PLUGIN=ON -DACTS_BUILD_TGEO_PLUGIN=ON -DACTS_BUILD_EXAMPLES=ON
  )
set_target_properties(external_acts PROPERTIES EXCLUDE_FROM_ALL TRUE)
message("To build ACTS package, use command:  make external_acts")
