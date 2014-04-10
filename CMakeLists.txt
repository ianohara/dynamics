cmake_minimum_required(VERSION 2.8)

project(mechanics)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Wall -Wextra -Wfatal-errors -Werror")

include_directories("${PROJECT_SOURCE_DIR}/inc")

add_subdirectory(tests)
add_subdirectory(src)
