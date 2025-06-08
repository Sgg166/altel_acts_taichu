
# string (REGEX MATCH "\\.el[1-9]" os_version_suffix ${CMAKE_SYSTEM})
# message("-- os_version_suffix:      ${os_version_suffix}")

cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(altel_acts VERSION 1.0.0 LANGUAGES CXX Fortran C)

message(STATUS "CMAKE_ROOT:                  ${CMAKE_ROOT}")
message(STATUS "CMAKE_VERSION:               ${CMAKE_VERSION}")
message(STATUS "CMAKE_HOST_SYSTEM:           ${CMAKE_HOST_SYSTEM}")
message(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "CMAKE_SYSTEM:                ${CMAKE_SYSTEM}")
message(STATUS "CMAKE_SYSTEM_PROCESSOR:      ${CMAKE_SYSTEM_PROCESSOR}")


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/INSTALL CACHE PATH "Default prefix of install directories" FORCE)
endif()
# set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME devel)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(CMAKE_LINK_SEARCH_START_STATIC ON)
# set(CMAKE_LINK_SEARCH_END_STATIC ON)

set(BUILD_SHARED_LIBS TRUE)

if(NOT BUILD_SHARED_LIBS)
  set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
endif()

set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_INSTALL_RPATH "@loader_path/../lib64;@loader_path/../lib")
else()
  set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib64:\$ORIGIN/../lib")
endif()
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# find_package order
set(CMAKE_FIND_USE_PACKAGE_ROOT_PATH       ON)   # 1. -DXX_ROOT and user ENV XX_ROOT
set(CMAKE_FIND_USE_CMAKE_PATH              ON)   # 2. -DCMAKE_PREFIX_PATH
set(CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH  OFF)  # 3. user ENV CMAKE_PREFIX_PATH  and XX_DIR
#                                                # 4. find_package option HINTS
set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH OFF)  # 5. ENV PATH
set(CMAKE_FIND_USE_PACKAGE_REGISTRY        OFF)  # 6. user registery: /home/user/.cmake
set(CMAKE_FIND_USE_CMAKE_SYSTEM_PATH       ON)   # 7. system CMAKE_SYSTEM_PREFIX_PATH
set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF)  # 8. system registery: somewhere
#                                                # 9. find_package option PATHS

## thread                                          
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# ## boost
if(EXISTS /usr/include/boost169)
  set(BOOST_INCLUDEDIR /usr/include/boost169)
endif()
if(EXISTS /usr/lib64/boost169)
  set(BOOST_LIBRARYDIR /usr/lib64/boost169)
endif()
set(Boost_NO_BOOST_CMAKE ON) # disable new cmake features from Boost 1.70 on
find_package(Boost 1.69 REQUIRED
  COMPONENTS
  filesystem program_options unit_test_framework
  # NO_CMAKE_SYSTEM_PATH
  )

## root
find_package(ROOT REQUIRED
  HINTS $ENV{ROOTSYS}/cmake
  # NO_CMAKE_SYSTEM_PATH
)


## external acts glfw
add_subdirectory(external)


if(ALTEL_EXTERNAL_CONFIG_ONLY)
  message("")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("!!!option ALTEL_EXTERNAL_CONFIG_ONLY is temporarily enabled,  (by -DALTEL_EXTERNAL_CONFIG_ONLY=ON)") 
  message("!!!cmake is now configured to build external packages only")
  message("!!!To build acts:")  
  message("!!!>> make external_acts")
  message("!!!To build eudaq:")  
  message("!!!>> make external_eudaq")
  message("!!!To build glfw:")  
  message("!!!>> make external_glfw")
  message("!!!Re-run cmake configuration will disable ALTEL_EXTERNAL_CONFIG_ONLY. Then you can continue to build the main software")
  message("!!!>> cmake ..")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("")
  message("")
  unset(ALTEL_EXTERNAL_CONFIG_ONLY CACHE)
  return()
endif()
  
find_package(glfw3 3.3
  HINTS
  ${CMAKE_INSTALL_PREFIX}/lib64/cmake/glfw3
)

if(NOT glfw3_FOUND)
  message("")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("!!!Missing glfw in system, please build glfw package by run command (AND then re-run cmake to complete config):")
  message("!!!>> make external_glfw")
  message("!!!>> cmake ..")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("")
  message("")
  return()
endif()

find_package(Acts
  COMPONENTS Core Alignment
  HINTS
  ${CMAKE_INSTALL_PREFIX}/cmake/Acts
  )

if(NOT Acts_FOUND)
  message("")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("!!!Missing Acts dependency, please build Acts package by run command (AND then re-run cmake to complete config): ")
  message("!!!>> make external_acts")
  message("!!!>> cmake ..")
  message("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  message("")
  message("")
  return()
endif()


get_filename_component(ACTS_INSTALL_DIR ${Acts_DIR}/../../.. ABSOLUTE)
set(ACTS_INSTALL_LIB_DIR  ${ACTS_INSTALL_DIR}/lib64)
message(STATUS "FOUND acts install cmake folder <Acts_DIR>:  ${Acts_DIR}")
message(STATUS "FOUND acts install folder <ACTS_INSTALL_DIR>:  ${ACTS_INSTALL_DIR}")

# include_directories(${ACTS_INSTALL_DIR}/include/ActsFatras)

add_subdirectory(source)

# install-tree
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../include")
configure_file(cmake/Config.cmake.in         ${PROJECT_BINARY_DIR}/cmake_install/${PROJECT_NAME}Config.cmake @ONLY)
configure_file(cmake/ConfigVersion.cmake.in  ${PROJECT_BINARY_DIR}/cmake_install/${PROJECT_NAME}ConfigVersion.cmake @ONLY)
install(FILES
  ${PROJECT_BINARY_DIR}/cmake_install/${PROJECT_NAME}Config.cmake
  ${PROJECT_BINARY_DIR}/cmake_install/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION cmake COMPONENT devel)

install(EXPORT ${PROJECT_NAME}Targets DESTINATION cmake)


# build-tree
set(CONF_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/../source/lib/include ${PROJECT_BINARY_DIR})
configure_file(cmake/Config.cmake.in          ${PROJECT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake @ONLY)
configure_file(cmake/ConfigVersion.cmake.in   ${PROJECT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake @ONLY)

export(PACKAGE ${PROJECT_NAME})

