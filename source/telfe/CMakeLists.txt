add_subdirectory(resource)


aux_source_directory(src THE_SRC_FILES)
add_library(altel-frontend SHARED ${THE_SRC_FILES} )
target_link_libraries (altel-frontend PUBLIC sensor-firmware-builtinJSON-headers mycommon stdc++fs altel-data)
target_include_directories(altel-frontend
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
)

set(THE_PUBLIC_HEADER
  include/Frontend.hh
  include/StreamInBuffer.hh
  include/TcpConnection.hh
  include/rbcp.hh
)

set_target_properties(altel-frontend PROPERTIES PUBLIC_HEADER "${THE_PUBLIC_HEADER}")

install(TARGETS altel-frontend
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
)

add_subdirectory(exe)

