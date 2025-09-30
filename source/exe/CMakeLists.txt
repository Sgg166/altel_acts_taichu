unset(EXE_TARGET_LIST)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore )

find_package(eudaq REQUIRED COMPONENTS eudaq::core)

set(ROOT_FEATURE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
get_target_property(ROOT_CXX_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_11" )
  set(ROOT_FEATURE_CXX_STANDARD 11)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_14" )
  set(ROOT_FEATURE_CXX_STANDARD 14)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_17" )
  set(ROOT_FEATURE_CXX_STANDARD 17)
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
     ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad
     altel-telfw altel-telgl galogen
     )
   


  add_executable(altelTelEventViewer altelTelEventViewer.cpp)
  list(APPEND EXE_TARGET_LIST altelTelEventViewer)
  target_include_directories(altelTelEventViewer  PRIVATE ./ )
  target_link_libraries(altelTelEventViewer
    altel-data-event altel-data-root
    mycommon
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  
    altel-telfw altel-telgl galogen
    )

  add_executable(altelEudaqRawViewer altelEudaqRawViewer.cpp)
  list(APPEND EXE_TARGET_LIST altelEudaqRawViewer)
  target_include_directories(altelEudaqRawViewer  PRIVATE ./ )
  target_link_libraries(altelEudaqRawViewer
    altel-data-event altel-data-eudaq
    mycommon
    altel-telfw altel-telgl galogen
    altel-data-root
    )


  add_executable(altelDataTaking altelDataTaking.cpp)
  list(APPEND EXE_TARGET_LIST altelDataTaking)
  target_include_directories(altelDataTaking  PRIVATE ./ ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(altelDataTaking
    altel-data-event altel-data-root altel-data-eudaq altel-rbcp
    ROOT::Core ROOT::RIO ROOT::Tree
    mycommon
    altel-telfw altel-telgl galogen
    )

endif()

add_executable(altelConvert altelConvert.cpp)
list(APPEND EXE_TARGET_LIST altelConvert)
target_link_libraries(altelConvert
  altel-acts altel-data-event altel-data-root
  mycommon
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  ROOT::MathCore
  )

add_executable(altelConvert2josn altelConvert2josn.cpp)
list(APPEND EXE_TARGET_LIST altelConvert2josn)
target_link_libraries(altelConvert2josn
  altel-acts altel-data-event 
  mycommon
  rapidjson
  eudaq::core
  )

add_executable(altelActsTrack altelActsTracking.cpp)
list(APPEND EXE_TARGET_LIST altelActsTrack)
target_link_libraries(altelActsTrack
  altel-acts altel-data-event altel-data-root
  mycommon
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  ROOT::MathCore
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
add_executable(test test.cc)
list(APPEND EXE_TARGET_LIST test)
target_include_directories(test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../teldata/event/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../telacts/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../teldaq/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../telfe/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
target_link_libraries(test
    PRIVATE
    altel-data-event altel-data-root altel-data-eudaq altel-rbcp
    altel-mille
    mycommon
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  ROOT::MathCore
    altel-telfw altel-telgl galogen
)


# add_executable(altelTcpClient altelTcpClient.cpp)
# list(APPEND EXE_TARGET_LIST altelTcpClient)
# target_link_libraries(altelTcpClient
#   PRIVATE
#   altel-frontend
#   altel-data-event
#   altel-data-root
#   altel-data-eudaq
#   altel-rbcp
#   mycommon
#   ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist  ROOT::MathCore
#   )


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
