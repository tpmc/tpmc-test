set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

add_executable(tpmc_test tpmc_test.cc)
target_link_libraries(tpmc_test ${TPMC_LIBRARIES})

add_executable(tpmc_test_1_tori tpmc_test_1_tori.cc)
target_link_libraries(tpmc_test_1_tori ${TPMC_LIBRARIES})
