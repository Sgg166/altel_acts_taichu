find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

set(ROOT_FEATURE_CXX_STANDARD ${CMAKE_CXX_STANDARD})
get_target_property(ROOT_CXX_COMPILE_FEATURES ROOT::Core INTERFACE_COMPILE_FEATURES)
if(ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_11" )
  set(ROOT_FEATURE_CXX_STANDARD 11)
endif()
if( ROOT_CXX_COMPILE_FEATURES STREQUAL "cxx_std_14" )
  set(ROOT_FEATURE_CXX_STANDARD 14)
endif()

aux_source_directory(src LIB_SRC)
add_library(altel-data-root SHARED ${LIB_SRC})

target_include_directories( altel-data-root
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  # $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )
target_link_libraries(
  altel-data-root
  PUBLIC altel-data-event
  ROOT::Core ROOT::RIO ROOT::Tree
  PRIVATE mycommon
  )

file(GLOB LIB_ALL_HEADERS "include/*.hpp")
set(LIB_PUBLIC_HEADERS "${LIB_ALL_HEADERS}" )
set_target_properties(altel-data-root PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

install(FILES macro/rootlogon.C  DESTINATION bin)

# install-tree
install(TARGETS altel-data-root
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

