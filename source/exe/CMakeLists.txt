
add_executable(AlignTelescopeGeometry AlignTelescopeGeometry.cpp)
target_link_libraries(AlignTelescopeGeometry PRIVATE
  altel-acts

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


add_executable(FitTelescopeTracks FitTelescopeTracks.cpp)
target_link_libraries(FitTelescopeTracks PRIVATE
  altel-acts

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

add_executable(EstimateTrackParameters TrackParametersEstimator.cpp)
target_link_libraries(EstimateTrackParameters PRIVATE
  altel-acts

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

# install-tree
install(TARGETS AlignTelescopeGeometry FitTelescopeTracks EstimateTrackParameters
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
