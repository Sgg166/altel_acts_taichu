
set(ACTS_SRC_DIR ${CMAKE_SOURCE_DIR}/../acts-core)
set(ACTS_INSTALL_DIR ${CMAKE_SOURCE_DIR}/../acts-core/INSTALL)

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



# set_target_properties(ActsExamplesDetectorsCommon PROPERTIES
#   IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesDetectorsCommon.so
#   INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Detectors/Common/include
#   )

# add_library(ActsExamplesPropagation SHARED IMPORTED)
# set_target_properties(ActsExamplesPropagation PROPERTIES
#   IMPORTED_LOCATION  ${ACTS_INSTALL_LIB_DIR}/libActsExamplesPropagation.so
#   INTERFACE_INCLUDE_DIRECTORIES ${ACTS_SRC_DIR}/Examples/Algorithms/Propagation/include
#   )

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
