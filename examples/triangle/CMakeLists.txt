add_executable(triangle main.cpp)

if(EXISTS ${BUILD_PATH}/lib${LIB_NAME}.a)
    target_link_libraries(triangle ${PROJECT_NAME})
endif()
