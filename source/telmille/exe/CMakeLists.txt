add_executable(TelMille_main TelMille_main.cpp)
target_link_libraries(TelMille_main PRIVATE mycommon altel-mille)

add_executable(GeometryFileUpdate_main GeometryFileUpdate_main.cpp)
target_link_libraries(GeometryFileUpdate_main PRIVATE mycommon altel-mille)

find_package(ROOT REQUIRED COMPONENTS RIO Net)

add_executable(LineFit_main LineFit_main.cpp)
target_link_libraries(LineFit_main PRIVATE mycommon altel-mille ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Graf3d)

install(TARGETS TelMille_main GeometryFileUpdate_main LineFit_main
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

