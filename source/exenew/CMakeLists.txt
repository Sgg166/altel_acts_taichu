
add_executable(TelDetectorAlign TelAlign.cpp)
list(APPEND EXE_TARGET_LIST TelDetectorAlign)
target_link_libraries(TelDetectorAlign
  altel-actsnew ActsAlignment mycommon)

add_executable(TelTrackFind TelTrackFind.cpp)
list(APPEND EXE_TARGET_LIST TelTrackFind)
target_link_libraries(TelTrackFind
  altel-actsnew ActsAlignment mycommon)

# install-tree
install(TARGETS ${EXE_TARGET_LIST}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
