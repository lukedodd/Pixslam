copy_files("copy_examples" ${CMAKE_CURRENT_SOURCE_DIR}/*.psm ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
   OUTPUT min_out.png
   COMMAND pixslam min.psm ../example_data/lena.png min_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT metaballs_out.png
   COMMAND pixslam metaballs.psm ../example_data/lena.png metaballs_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT metaballs_binary_out.png
   COMMAND pixslam metaballs_binary.psm ../example_data/lena.png metaballs_binary_out.png
   DEPENDS pixslam
)

ADD_CUSTOM_TARGET(examples ALL DEPENDS copy_examples copy_example_data min_out.png metaballs_out.png metaballs_binary_out.png)

