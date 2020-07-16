cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(altel_acts VERSION 1.0.0 LANGUAGES CXX)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/INSTALL CACHE PATH "Default prefix of install directories" FORCE)
endif()
# set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME devel)

set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# set(CMAKE_LINK_SEARCH_START_STATIC ON)
# set(CMAKE_LINK_SEARCH_END_STATIC ON)

# set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")


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


                                                  
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_USE_STATIC_LIBS ON) 
find_package(Boost REQUIRED COMPONENTS system thread filesystem chrono regex program_options unit_test_framework)

find_package(ROOT REQUIRED HINTS $ENV{ROOTSYS}/cmake NO_CMAKE_SYSTEM_PATH)


find_package(Acts REQUIRED COMPONENTS  Core Fatras DigitizationPlugin IdentificationPlugin

  HINTS
  ${PROJECT_SOURCE_DIR}/../acts-core/INSTALL/share/cmake/Acts
  )
# find_package(Acts REQUIRED COMPONENTS  JsonPlugin IdentificationPluginTGeoPlugin)


include(cmake/ActsExamples.cmake)
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

