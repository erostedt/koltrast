cmake_minimum_required(VERSION 3.25)

project(koltrast)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS
    "-Wall -pedantic -Wextra -Werror -Wshadow -Wsign-conversion -Wconversion -Wunreachable-code -g"
)

include_directories(src)
add_executable(${PROJECT_NAME} src/main.cpp)
