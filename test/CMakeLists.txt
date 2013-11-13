set(PROJECT_TEST_NAME "fit_gtr_test")
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(${PROJECT_TEST_NAME} "example_test.cpp")
target_link_libraries(${PROJECT_TEST_NAME} gtest gtest_main)

add_test(${PROJECT_TEST_NAME} ${PROJECT_TEST_NAME})
