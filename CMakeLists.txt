cmake_minimum_required(VERSION 2.8)

project(mechanics)

include_directories("${PROJECT_SOURCE_DIR}/inc")

add_subdirectory(tests)
add_subdirectory(src)
