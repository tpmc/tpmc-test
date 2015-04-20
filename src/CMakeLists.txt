set(CMAKE_BUILD_TYPE Release)

add_executable(tpmc_test tpmc_test.cc)
target_link_libraries(tpmc_test ${TPMC_LIBRARIES})

add_executable(tpmc_test_1_tori tpmc_test_1_tori.cc)
target_link_libraries(tpmc_test_1_tori ${TPMC_LIBRARIES})

add_executable(tpmc_test_2_rotation tpmc_test_2_rotation.cc)
target_link_libraries(tpmc_test_2_rotation ${TPMC_LIBRARIES})

add_executable(tpmc_test_3_performance tpmc_test_3_performance.cc)
target_link_libraries(tpmc_test_3_performance ${TPMC_LIBRARIES})
