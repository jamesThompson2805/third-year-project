cmake_minimum_required(VERSION 3.26)

project(third_year_project)
enable_testing()


set(CMAKE_CXX_VERSION 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(${PROJECT_NAME} src/main.cpp
			       src/plotting/series_plotting.cpp
			       src/plotting/plot_types.cpp
			       src/evaluations/general.cpp
			       src/evaluations/capla.cpp)

add_subdirectory(lib/sequence_gen)
add_subdirectory(lib/parsing)
add_subdirectory(lib/dimension_reductions)
add_subdirectory(lib/similarity_search)
add_subdirectory(lib/cleaning)
target_link_libraries(${PROJECT_NAME} PUBLIC my_sequence_gen)
target_link_libraries(${PROJECT_NAME} PUBLIC my_parsing)
target_link_libraries(${PROJECT_NAME} PUBLIC my_dimension_reductions)
target_link_libraries(${PROJECT_NAME} PUBLIC my_similarity_search)
target_link_libraries(${PROJECT_NAME} PUBLIC my_cleaning)

add_subdirectory(external/libs/gnuplot-iostream)
target_link_libraries(${PROJECT_NAME} PUBLIC gnuplot_iostream)
add_subdirectory(external/libs/pgbar)
target_link_libraries(${PROJECT_NAME} PUBLIC konvt_pgbar)


find_package(GTest REQUIRED)
set(TEST_NAME "third_year_testing")
add_executable(${TEST_NAME}
  tst/similarity_search/r-tree-test.cpp
  tst/dimension_reductions/dac_curve_fitting_test.cpp
  tst/dimension_reductions/double_window_test.cpp)
target_link_libraries( ${TEST_NAME} PUBLIC GTest::gtest_main)
target_link_libraries(${TEST_NAME} PUBLIC my_sequence_gen)
target_link_libraries(${TEST_NAME} PUBLIC my_parsing)
target_link_libraries(${TEST_NAME} PUBLIC my_dimension_reductions)
target_link_libraries(${TEST_NAME} PUBLIC my_similarity_search)
target_link_libraries(${TEST_NAME} PUBLIC my_cleaning)

include(GoogleTest)
gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)
add_test(${TEST_NAME} ${TEST_NAME})
