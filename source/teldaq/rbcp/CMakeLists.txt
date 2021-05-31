add_subdirectory(rbcp)

aux_source_directory(src LIB_SRC)
add_library(altel-rbcp STATIC ${LIB_SRC} )
target_link_libraries(altel-rbcp PUBLIC mycommon  altel-data-event PRIVATE rbcp)

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

create_cpp_raw_string(
  ${CMAKE_CURRENT_BINARY_DIR}/include/altel_reg_cmd_list.hh
  ${CMAKE_CURRENT_LIST_DIR}/misc/altel_reg_cmd_list.json
  )


target_include_directories(altel-rbcp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  )

install(TARGETS altel-rbcp
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib      COMPONENT runtime
  ARCHIVE       DESTINATION lib      COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
