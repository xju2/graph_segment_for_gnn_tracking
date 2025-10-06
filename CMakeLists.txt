cmake_minimum_required(VERSION 3.16...3.31)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

include(FetchContent)

# Option to force vendored Boost (default OFF for dev, ON for wheels)
option(GRAPH_SEGMENT_VENDOR_BOOST "Use vendored Boost instead of system Boost"
       OFF)

if(NOT GRAPH_SEGMENT_VENDOR_BOOST)
  find_package(Boost QUIET COMPONENTS graph)
endif()

if(NOT Boost_FOUND)
  message(STATUS "Fetching Boost from source...")

  FetchContent_Declare(
    boost
    GIT_REPOSITORY https://github.com/boostorg/boost.git
    GIT_TAG boost-1.89.0)

  set(BOOST_ENABLE_CMAKE
      ON
      CACHE BOOL "Enable CMake in Boost")
  FetchContent_MakeAvailable(boost)

  set(Boost_INCLUDE_DIRS ${boost_SOURCE_DIR})

endif()

message(STATUS "Boost Libraries: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")

add_subdirectory(src)
