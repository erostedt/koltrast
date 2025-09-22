cmake_minimum_required(VERSION 3.25)

project(koltrast)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS
    "-Wall -pedantic -Wextra -Werror -Wshadow -Wsign-conversion -Wconversion -Wunreachable-code -g"
)

set(EXAMPLES_PATH examples)
set(BUILD_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build)
set(LIB_NAME ${PROJECT_NAME})
set(SRC_PATH src)
set(THIRD_PARTY third_party)

include_directories(${SRC_PATH})
file(GLOB source_files ${SRC_PATH}/*.cpp)

add_library(${LIB_NAME} ${source_files})
target_include_directories(${LIB_NAME} SYSTEM PUBLIC ${THIRD_PARTY})

target_link_libraries(${LIB_NAME} tbb)

add_subdirectory(tests)
add_subdirectory(examples)
