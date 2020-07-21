
set(Boost_NO_BOOST_CMAKE ON) # disable new cmake features from Boost 1.70 on

if(EXISTS /usr/include/boost169)
  set(BOOST_INCLUDEDIR /usr/include/boost169)
endif()
  
find_package(Boost 1.69 REQUIRED
  COMPONENTS
  system thread filesystem chrono regex program_options unit_test_framework
  # NO_CMAKE_SYSTEM_PATH
  )

find_package(ROOT REQUIRED HINTS $ENV{ROOTSYS}/cmake NO_CMAKE_SYSTEM_PATH)

if(NOT DEFINED ACTS_SRC_DIR)
  set(ACTS_SRC_DIR ACTS_SRC_DIR-NOTFOUND  CACHE PATH "ACTS top source folder")
  message(FATAL_ERROR "<ACTS_SRC_DIR> is not defined, please set it to ACTS top source folder")
endif()

find_package(Acts REQUIRED COMPONENTS  Core Fatras DigitizationPlugin IdentificationPlugin TGeoPlugin # JsonPlugin
  HINTS
  ${ACTS_SRC_DIR}/INSTALL/share/cmake/Acts
  )

set(ACTS_INSTALL_DIR ${Acts_DIR}/../../..)

set(ACTS_INSTALL_LIB_DIR  ${ACTS_INSTALL_DIR}/lib64)

include_directories(${ACTS_INSTALL_DIR}/include/ActsFatras)

add_library(ActsExamplesFramework SHARED IMPORTED)
set_target_properties(ActsExamplesFramework PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesFramework.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Framework/include
  )

add_library(ActsExamplesCommon SHARED IMPORTED)
set_target_properties(ActsExamplesCommon PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesCommon.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Run/Common/include
  )

add_library(ActsExamplesDetectorsCommon INTERFACE)
target_include_directories(ActsExamplesDetectorsCommon INTERFACE
  ${ACTS_SRC_DIR}/Examples/Detectors/Common/include
  ${ACTS_SRC_DIR}/Examples/Algorithms/Propagation/include
  )


add_library(ActsExamplesPropagation INTERFACE)
target_include_directories(ActsExamplesPropagation INTERFACE
  ${ACTS_SRC_DIR}/Examples/Algorithms/Propagation/include
  )

add_library(ActsExamplesTruthTracking SHARED IMPORTED)
set_target_properties(ActsExamplesTruthTracking PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesTruthTracking.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Algorithms/TruthTracking
  )

add_library(ActsExamplesMagneticField SHARED IMPORTED)
set_target_properties(ActsExamplesMagneticField PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesMagneticField.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Detectors/MagneticField/include
  )

add_library(ActsExamplesAlignment SHARED IMPORTED)
set_target_properties(ActsExamplesAlignment PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesAlignment.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Algorithms/Alignment/include
  )

add_library(ActsExamplesIoPerformance SHARED IMPORTED)
set_target_properties(ActsExamplesIoPerformance PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesIoPerformance.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Io/Performance
  )

add_library(ActsExamplesIoRoot SHARED IMPORTED)
set_target_properties(ActsExamplesIoRoot PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesIoRoot.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Io/Root
  )

add_library(ActsExamplesIoObj SHARED IMPORTED)
set_target_properties(ActsExamplesIoObj PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesIoObj.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Io/Obj
  )

add_library(ActsExamplesIoCsv SHARED IMPORTED)
set_target_properties(ActsExamplesIoCsv PROPERTIES
  IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesIoCsv.so
  INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Io/Csv
  )
