add_executable(x11 main.cpp window.hpp)

target_link_libraries(x11 ${PROJECT_NAME} X11)
