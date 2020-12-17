
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

set(ROOT_FEATURE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
get_target_property(ROOT_CXX_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_11" )
  set(ROOT_FEATURE_CXX_STANDARD 11)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_14" )
  set(ROOT_FEATURE_CXX_STANDARD 14)
endif()


if(TARGET altel-telfw)
  add_executable(TelDetectorResidual TelDetectorResidual.cpp)
  list(APPEND EXE_TARGET_LIST TelDetectorResidual)
  target_include_directories(TelDetectorResidual  PRIVATE ./ )
  target_link_libraries(TelDetectorResidual
    altel-acts ActsAlignment mycommon
    ROOT::Core ROOT::RIO ROOT::Tree
    altel-telfw altel-telgl galogen
    )


  add_executable(TelKinkAngle TelKinkAngle.cpp)
  list(APPEND EXE_TARGET_LIST TelKinkAngle)
  target_include_directories(TelKinkAngle  PRIVATE ./ )
  target_link_libraries(TelKinkAngle
    altel-acts ActsAlignment mycommon
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
    altel-telfw altel-telgl galogen
    )

  add_executable(altelTelEventViewer altelTelEventViewer.cpp)
  list(APPEND EXE_TARGET_LIST altelTelEventViewer)
  target_include_directories(altelTelEventViewer  PRIVATE ./ )
  target_link_libraries(altelTelEventViewer
    altel-data-event altel-data-root
    mycommon
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  
    )

endif()


add_executable(altelActsTrack altelActsTracking.cpp)
list(APPEND EXE_TARGET_LIST altelActsTrack)
target_link_libraries(altelActsTrack
  altel-acts altel-data-event altel-data-root
  mycommon
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  
  )

add_executable(altelMilleBin altelMilleAlign.cpp)
list(APPEND EXE_TARGET_LIST altelMilleBin)
target_include_directories(altelMilleBin  PRIVATE ./ )
target_link_libraries(altelMilleBin
  PRIVATE
  altel-data-event
  altel-data-eudaq
  altel-mille
  mycommon
  )



install(FILES altelTTreeExampleReader.C DESTINATION bin)

# install-tree
install(TARGETS ${EXE_TARGET_LIST}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
