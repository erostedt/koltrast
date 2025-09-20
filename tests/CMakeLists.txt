file(GLOB test_files *.cpp *hpp)

add_executable(run_tests ${test_files})

target_link_libraries(run_tests ${PROJECT_NAME})
