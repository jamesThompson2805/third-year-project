cmake_minimum_required(VERSION 3.26)

project(my_dimension_reductions)

set(CMAKE_CXX_STANDARD 17)

add_library( ${PROJECT_NAME} naive_pla.cpp paa.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})



