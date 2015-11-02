project(tpmc_test)
cmake_minimum_required(VERSION 2.8.6)

enable_testing()

set(TPMC_UPSTREAM "https://www.github.com/tpmc/tpmc.git" CACHE STRING "git URL of the TPMC library")
set(TPMC_TEST_DIR "${CMAKE_SOURCE_DIR}/src/full_tests/" CACHE STRING "path to the test configuration")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_TPMC_PROFILING -DNDEBUG -std=c++11 -O3 -march=native -funroll-loops -pedantic -Wall -Wextra -DTPMC_TEST_DIR=\"${TPMC_TEST_DIR}\"")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(external)
include_directories(${TPMC_INCLUDE_DIR})
add_subdirectory(src)
