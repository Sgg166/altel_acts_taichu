
add_executable(TelDetectorAlign TelAlign.cpp)
list(APPEND EXE_TARGET_LIST TelDetectorAlign)
target_link_libraries(TelDetectorAlign
  altel-acts ActsAlignment mycommon)

add_executable(TelTrackFind TelTrackFind.cpp)
list(APPEND EXE_TARGET_LIST TelTrackFind)
target_link_libraries(TelTrackFind
  altel-acts ActsAlignment mycommon)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

set(ROOT_FEATURE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
get_target_property(ROOT_CXX_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_11" )
  set(ROOT_FEATURE_CXX_STANDARD 11)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_14" )
  set(ROOT_FEATURE_CXX_STANDARD 14)
endif()


add_executable(TelTrackResidual TelTrackResidual.cpp)
list(APPEND EXE_TARGET_LIST TelTrackResidual)
target_link_libraries(TelTrackResidual
  altel-acts ActsAlignment mycommon
  ROOT::Core ROOT::RIO ROOT::Tree
  )


add_executable(TelDetectorResidual TelDetectorResidual.cpp)
list(APPEND EXE_TARGET_LIST TelDetectorResidual)
target_include_directories(TelDetectorResidual  PRIVATE ./ )
target_link_libraries(TelDetectorResidual
  altel-acts ActsAlignment mycommon
  ROOT::Core ROOT::RIO ROOT::Tree
  altel-telfw altel-telgl galogen
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
