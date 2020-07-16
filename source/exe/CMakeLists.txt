
add_executable(RecTelescopeTracks RecTelescopeTracks.cpp)
target_link_libraries(RecTelescopeTracks PRIVATE
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
install(TARGETS RecTelescopeTracks
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
