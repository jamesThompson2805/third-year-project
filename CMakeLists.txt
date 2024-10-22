cmake_minimum_required(VERSION 3.26)

project(third_year_project)

set(CMAKE_CXX_VERSION 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(${PROJECT_NAME} src/main.cpp)

add_subdirectory(lib/sequence_gen)
add_subdirectory(lib/parsing)
add_subdirectory(lib/dimension_reductions)
add_subdirectory(lib/similarity_search)
target_link_libraries(${PROJECT_NAME} PUBLIC my_sequence_gen)
target_link_libraries(${PROJECT_NAME} PUBLIC my_parsing)
target_link_libraries(${PROJECT_NAME} PUBLIC my_dimension_reductions)
target_link_libraries(${PROJECT_NAME} PUBLIC my_similarity_search)
