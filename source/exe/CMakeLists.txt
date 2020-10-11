
add_executable(AlignTelescopeGeometry AlignTelescopeGeometry.cpp)
list(APPEND EXE_TARGET_LIST AlignTelescopeGeometry)
target_link_libraries(AlignTelescopeGeometry PRIVATE
  altel-acts

  ActsAlignment

  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesDetectorsCommon
  ActsExamplesMagneticField
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoPerformance
  rapidjson

  ROOT::Core
  ROOT::Hist
  ROOT::Tree
  )


add_executable(EstimateTrackParameters EstimateTrackParameters.cpp)
list(APPEND EXE_TARGET_LIST EstimateTrackParameters)
target_link_libraries(EstimateTrackParameters PRIVATE
  altel-acts

  ActsAlignment

  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesDetectorsCommon
  ActsExamplesMagneticField
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoPerformance
  rapidjson

  ROOT::Core
  ROOT::Hist
  ROOT::Tree
  Boost::boost
  )



add_executable(FindTelescopeTracks FindTelescopeTracks.cpp)
list(APPEND EXE_TARGET_LIST FindTelescopeTracks)
target_link_libraries(FindTelescopeTracks PRIVATE
  altel-acts

  ActsAlignment

  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesDetectorsCommon
  ActsExamplesMagneticField
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoPerformance
  rapidjson

  ROOT::Core
  ROOT::Hist
  ROOT::Tree
  )


add_executable(PrintEventJson PrintEventJson.cpp)
list(APPEND EXE_TARGET_LIST PrintEventJson)
target_link_libraries(PrintEventJson PRIVATE
  altel-acts
  rapidjson
  )




# install-tree
install(TARGETS ${EXE_TARGET_LIST}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
