cmake_minimum_required(VERSION 3.26)

project(my_similarity_search)

set(CMAKE_CXX_STANDARD 17)

add_library( ${PROJECT_NAME} sequential_scan.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})



