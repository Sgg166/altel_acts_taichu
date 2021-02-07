
unset(EXE_TARGET_LIST)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree MathCore)


function(create_cpp_raw_string)
  if( (${ARGC} LESS 2) OR (${ARGC} GREATER 3 ) )
    message(FATAL_ERROR "create_cpp_raw_string require 2 or 3 arguments, but given ${ARGC}")
  endif()
  if(${ARGC} EQUAL 2)
    set(RAW_STRING_DELIM "~~~cpp~raw~~~")
  endif()
  if(${ARGC} EQUAL 3)
    set(RAW_STRING_DELIM ${ARGV2})
  endif()
  if( (NOT EXISTS ${ARGV0}) OR (${ARGV1} IS_NEWER_THAN ${ARGV0}) )
    file(READ ${ARGV1} RAW_STRING_CONTENT)
    set(RAW_STRING_CONTENT "R\"${RAW_STRING_DELIM}(${RAW_STRING_CONTENT})${RAW_STRING_DELIM}\"")
    file(WRITE ${ARGV0} "${RAW_STRING_CONTENT}")
  endif()
endfunction()

set(ROOT_FEATURE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
get_target_property(ROOT_CXX_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_11" )
  set(ROOT_FEATURE_CXX_STANDARD 11)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_14" )
  set(ROOT_FEATURE_CXX_STANDARD 14)
endif()

if(TARGET altel-telfw)
  # add_executable(TelDetectorResidual TelDetectorResidual.cpp)
  # list(APPEND EXE_TARGET_LIST TelDetectorResidual)
  # target_include_directories(TelDetectorResidual  PRIVATE ./ )
  # target_link_libraries(TelDetectorResidual
  #   altel-acts ActsAlignment mycommon
  #   ROOT::Core ROOT::RIO ROOT::Tree
  #   altel-telfw altel-telgl galogen
  #   )

  # add_executable(TelKinkAngle TelKinkAngle.cpp)
  # list(APPEND EXE_TARGET_LIST TelKinkAngle)
  # target_include_directories(TelKinkAngle  PRIVATE ./ )
  # target_link_libraries(TelKinkAngle
  #   altel-acts ActsAlignment mycommon
  #   ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
  #   altel-telfw altel-telgl galogen
  #   )

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
    )

  create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/altel_run_conf.hh
    ${PROJECT_SOURCE_DIR}/source/teldaq/rbcp/misc/altel_run_conf.json)

  add_executable(altelDataTaking altelDataTaking.cpp)
  list(APPEND EXE_TARGET_LIST altelDataTaking)
  target_include_directories(altelDataTaking  PRIVATE ./ ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(altelDataTaking
    altel-data-event altel-data-eudaq altel-rbcp
    mycommon
    altel-telfw altel-telgl galogen
    )

endif()


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
