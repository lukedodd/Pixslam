cmake_minimum_required(VERSION 2.8)
project(PIXSLAM)

add_subdirectory(libs)
add_subdirectory(examples)

# C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 

include_directories(${PIXSLAM_SOURCE_DIR}/src ${PIXSLAM_SOURCE_DIR}/libs/asmjit/src ${PIXSLAM_SOURCE_DIR}/libs/stb_image)

# aux_source_directory(src ${SOURCES})
add_executable(pixslam src/main.cpp)

target_link_libraries(pixslam asmjit stb_image)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/example_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
