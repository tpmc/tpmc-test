project(tpmc_test)
cmake_minimum_required(VERSION 2.8.6)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
list(APPEND CMAKE_CXX_FLAGS "-std=c++11")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Tpmc REQUIRED)
include_directories(${TPMC_INCLUDE_DIRS})

add_subdirectory(src)
