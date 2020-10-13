
add_executable(TelAlign TelAlign.cpp)
list(APPEND EXE_TARGET_LIST TelAlign)
target_link_libraries(TelAlign
  altel-actsnew
  ActsAlignment
  mycommon
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
