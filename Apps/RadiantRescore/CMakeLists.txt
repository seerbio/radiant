cmake_minimum_required(VERSION 3.18)
project(RadiantRescore)
find_package(Qt5 COMPONENTS Core REQUIRED)
add_executable(RadiantRescore main.cpp)
target_link_libraries(RadiantRescore Qt5::Core AlgorithmsFFLib FileReadersLib)
set_target_properties(RadiantRescore PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    CXX_STANDARD 17 CXX_EXTENSIONS OFF)
