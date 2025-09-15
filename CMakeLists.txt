cmake_minimum_required(VERSION 3.25)

project(koltrast)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS
    "-Wall -pedantic -Wextra -Werror -Wshadow -Wsign-conversion -Wconversion -Wunreachable-code -g"
)

set(EXAMPLES_PATH examples)
set(BUILD_PATH build)
set(LIB_NAME ${PROJECT_NAME})
set(SRC_PATH src)
set(THIRD_PARTY third_party)

include_directories(${SRC_PATH})
include_directories(${THIRD_PARTY})
file(GLOB source_files ${SRC_PATH}/*.cpp)

if(source_files)
  add_library(${LIB_NAME} ${source_files})
endif()

add_subdirectory(tests)
add_subdirectory(examples)
