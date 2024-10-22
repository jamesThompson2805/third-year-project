cmake_minimum_required(VERSION 3.26)

project(my_parsing)

set(CMAKE_CXX_STANDARD 17)

add_library( ${PROJECT_NAME} tsv_parsing.cpp ucr_parsing.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})



