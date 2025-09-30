# - Config file for the package
# It defines the following variables
#  altel_acts_INCLUDE_DIRS - include directories
#  altel_acts_LIBRARIES    - libraries to link against
#  altel_acts_EXECUTABLE   - the  executable

# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET altel-acts)
  include(${CMAKE_CURRENT_LIST_DIR}/altel_actsTargets.cmake)
endif()

# These are IMPORTED targets created by FooBarTargets.cmake
# set("altel_acts"_LIBRARIES )
# set("altel_acts"_EXECUTABLE )

set(altel_acts_LIBRARIES altel-acts)
set(altel_acts_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../include")
