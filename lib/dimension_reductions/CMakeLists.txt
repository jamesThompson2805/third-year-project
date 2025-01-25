cmake_minimum_required(VERSION 3.26)

project(my_dimension_reductions)

set(CMAKE_CXX_STANDARD 17)
enable_testing()

add_library( ${PROJECT_NAME} pla.cpp paa.cpp dac_curve_fitting.cpp double_window.cpp exact_dp.cpp conv_double_window.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})
