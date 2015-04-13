set(CMAKE_BUILD_TYPE Release)
add_executable(tpmc_test tpmc_test.cc)
target_link_libraries(tpmc_test ${TPMC_LIBRARIES})
