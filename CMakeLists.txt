cmake_minimum_required(VERSION 3.25)

project(koltrast)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(EXAMPLES_PATH examples)
set(LIB_PATH trend)
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_FLAGS
    "-Wall -pedantic -Wextra -Werror -Wshadow -Wsign-conversion -Wconversion -Wunreachable-code -g"
)

add_executable(${PROJECT_NAME} main.cpp)
