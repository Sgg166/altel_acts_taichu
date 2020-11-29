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

create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/TelVertex_glsl.hh TelVertex.glsl)
create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/TelGeometry_glsl.hh TelGeometry.glsl)
create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/TelFragment_glsl.hh TelFragment.glsl)
create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/HitVertex_glsl.hh HitVertex.glsl)
create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/HitGeometry_glsl.hh HitGeometry.glsl)
create_cpp_raw_string(${CMAKE_CURRENT_BINARY_DIR}/include/TrackVertex_glsl.hh TrackVertex.glsl)

add_library(altel-glsl INTERFACE)
target_include_directories(altel-glsl
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/include
  )
