file(GLOB test_files *.cpp *hpp)

add_executable(run_tests ${test_files})

if(EXISTS ${BUILD_PATH}/lib${LIB_NAME}.a)
  target_link_libraries(run_tests ${PROJECT_NAME})
endif()
