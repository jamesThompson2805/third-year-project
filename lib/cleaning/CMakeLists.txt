cmake_minimum_required(VERSION 3.26)

project(my_cleaning)

set(CMAKE_CXX_STANDARD 17)

add_library( ${PROJECT_NAME} z_norm.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})



