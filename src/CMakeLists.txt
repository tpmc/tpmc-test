add_library(tpmc_test timer.cc geometry.cc)

add_executable(tpmc_test_1_tori tpmc_test_1_tori.cc)
target_link_libraries(tpmc_test_1_tori ${TPMC_LIBRARIES})
target_link_libraries(tpmc_test_1_tori tpmc_test)
add_dependencies(tpmc_test tpmc)

add_executable(tpmc_test_2_rotation tpmc_test_2_rotation.cc)
target_link_libraries(tpmc_test_2_rotation ${TPMC_LIBRARIES})
target_link_libraries(tpmc_test_2_rotation tpmc_test)
add_dependencies(tpmc_test tpmc)

add_executable(tpmc_test_3_performance tpmc_test_3_performance.cc)
target_link_libraries(tpmc_test_3_performance ${TPMC_LIBRARIES})
target_link_libraries(tpmc_test_3_performance tpmc_test)
add_dependencies(tpmc_test tpmc)

add_test(test1 tpmc_test_1_tori basic)
add_test(test2 tpmc_test_2_rotation basic)
add_test(test3 tpmc_test_3_performance basic)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tpmc_test_1_tori tpmc_test_2_rotation tpmc_test_3_performance)
