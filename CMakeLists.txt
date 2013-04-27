cmake_minimum_required(VERSION 2.8)
project(PIXSLAM)

macro(copy_files TARGET_NAME GLOBPAT DESTINATION)
  file(GLOB COPY_FILES
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    ${GLOBPAT})
  add_custom_target(${TARGET_NAME} ALL
    COMMENT "Copying files: ${GLOBPAT}")

  foreach(FILENAME ${COPY_FILES})
    set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
    set(DST "${DESTINATION}/${FILENAME}")

    add_custom_command(
      TARGET ${TARGET_NAME}
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
      )
  endforeach(FILENAME)
endmacro(copy_files)

add_subdirectory(libs)
add_subdirectory(examples)
add_subdirectory(example_data)

# C++11
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 

include_directories(${PIXSLAM_SOURCE_DIR}/src ${PIXSLAM_SOURCE_DIR}/libs/asmjit/src ${PIXSLAM_SOURCE_DIR}/libs/stb_image)

# aux_source_directory(src ${SOURCES})
add_executable(pixslam src/main.cpp)

target_link_libraries(pixslam asmjit stb_image)
