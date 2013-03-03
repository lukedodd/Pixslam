cmake_minimum_required(VERSION 2.8)
project(PIXSLAM)

add_subdirectory(libs)

# C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 

include_directories(${PIXSLAM_SOURCE_DIR}/src ${PIXSLAM_SOURCE_DIR}/libs/asmjit/src)

# aux_source_directory(src ${SOURCES})
add_executable(pixslam src/main.cpp)

target_link_libraries(pixslam asmjit stb_image)
