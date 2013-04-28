copy_files("copy_examples" ${CMAKE_CURRENT_SOURCE_DIR}/*.psm ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
   OUTPUT compose_out.png
   COMMAND pixslam compose.psm ../example_data/lena.png ../example_data/duck.png compose_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT erode_3x3_out.png
   COMMAND pixslam erode_3x3.psm ../example_data/lena.png erode_3x3_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT dilate_3x3_out.png
   COMMAND pixslam dilate_3x3.psm ../example_data/lena.png dilate_3x3_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT min_out.png
   COMMAND pixslam min.psm ../example_data/lena.png ../example_data/duck.png min_out.png
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

add_custom_command(
   OUTPUT metaballs_fancy_out.png
   COMMAND pixslam metaballs_fancy.psm ../example_data/lena.png metaballs_fancy_out.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT box_3x3_out_1.png
   COMMAND pixslam box_3x3.psm ../example_data/lena.png box_3x3_out_1.png
   DEPENDS pixslam
)

add_custom_command(
   OUTPUT box_3x3_out_2.png
   COMMAND pixslam box_3x3.psm metaballs_fancy_out.png box_3x3_out_2.png
   DEPENDS pixslam
)

ADD_CUSTOM_TARGET(examples ALL DEPENDS copy_examples copy_example_data compose_out.png erode_3x3_out.png dilate_3x3_out.png min_out.png metaballs_out.png metaballs_binary_out.png metaballs_fancy_out.png box_3x3_out_2.png)

