copy_files("copy_examples" ${CMAKE_CURRENT_SOURCE_DIR}/*.psm ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
   OUTPUT compose_out.png
   COMMAND pixslam --logCommand compose.psm ../example_data/lena.png ../example_data/duck.png compose_out.png
   DEPENDS pixslam compose.psm
)

add_custom_command(
   OUTPUT erode_3x3_out.png
   COMMAND pixslam --logCommand erode_3x3.psm ../example_data/lena.png erode_3x3_out.png
   DEPENDS pixslam erode_3x3.psm
)

add_custom_command(
   OUTPUT dilate_3x3_out.png
   COMMAND pixslam --logCommand dilate_3x3.psm ../example_data/lena.png dilate_3x3_out.png
   DEPENDS pixslam dilate_3x3.psm
)

add_custom_command(
   OUTPUT erode_5x5_out.png
   COMMAND pixslam --logCommand erode_5x5.psm ../example_data/lena.png erode_5x5_out.png
   DEPENDS pixslam erode_5x5.psm
)

add_custom_command(
   OUTPUT dilate_5x5_out.png
   COMMAND pixslam --logCommand dilate_5x5.psm ../example_data/lena.png dilate_5x5_out.png
   DEPENDS pixslam dilate_5x5.psm
)

add_custom_command(
   OUTPUT min_out.png
   COMMAND pixslam --logCommand min.psm ../example_data/lena.png ../example_data/duck.png min_out.png
   DEPENDS pixslam min.psm
)

add_custom_command(
   OUTPUT metaballs_out.png
   COMMAND pixslam --logCommand metaballs.psm ../example_data/lena.png metaballs_out.png
   DEPENDS pixslam metaballs.psm
)

add_custom_command(
   OUTPUT metaballs_binary_out.png
   COMMAND pixslam --logCommand metaballs_binary.psm ../example_data/lena.png metaballs_binary_out.png
   DEPENDS pixslam metaballs_binary.psm
)

add_custom_command(
   OUTPUT metaballs_fancy_out.png
   COMMAND pixslam --logCommand metaballs_fancy.psm ../example_data/lena.png metaballs_fancy_out.png
   DEPENDS pixslam metaballs_fancy.psm
)

add_custom_command(
   OUTPUT box_3x3_out_1.png
   COMMAND pixslam --logCommand box_3x3.psm ../example_data/lena.png box_3x3_out_1.png
   DEPENDS pixslam box_3x3.psm
)

add_custom_command(
   OUTPUT box_3x3_out_2.png
   COMMAND pixslam --logCommand box_3x3.psm metaballs_fancy_out.png box_3x3_out_2.png
   DEPENDS pixslam box_3x3.psm metaballs_fancy_out.png
)

add_custom_command(
   OUTPUT box_5x5_out_1.png
   COMMAND pixslam --logCommand box_5x5.psm ../example_data/lena.png box_5x5_out_1.png
   DEPENDS pixslam box_5x5.psm
)

add_custom_command(
   OUTPUT box_5x5_out_2.png
   COMMAND pixslam --logCommand box_5x5.psm metaballs_fancy_out.png box_5x5_out_2.png
   DEPENDS pixslam box_5x5.psm metaballs_fancy_out.png
)

add_custom_command(
   OUTPUT lena_edge.png
   COMMAND pixslam --logCommand absdiff.psm ../example_data/lena.png box_5x5_out_1.png lena_edge.png
   DEPENDS pixslam absdiff.psm box_5x5_out_2.png 
)

add_custom_command(
   OUTPUT metaball_edges.png
   COMMAND pixslam --logCommand absdiff.psm metaballs_fancy_out.png box_5x5_out_2.png metaball_edges.png
   DEPENDS pixslam absdiff.psm metaballs_fancy_out.png box_5x5_out_2.png 
)

add_custom_command(
   OUTPUT flipped_duck.png
   COMMAND pixslam --logCommand flip_vertical.psm ../example_data/duck.png flipped_duck.png
   DEPENDS pixslam flip_vertical.psm 
)

add_custom_command(
   OUTPUT threshold.png
   COMMAND pixslam --logCommand threshold.psm ../example_data/lena.png threshold.png
   DEPENDS pixslam threshold.psm 
)

add_custom_target(examples ALL DEPENDS copy_examples copy_example_data compose_out.png erode_3x3_out.png dilate_3x3_out.png erode_5x5_out.png dilate_5x5_out.png min_out.png metaballs_out.png metaballs_binary_out.png metaballs_fancy_out.png box_3x3_out_1.png box_3x3_out_2.png box_5x5_out_1.png box_5x5_out_2.png metaball_edges.png lena_edge.png flipped_duck.png threshold.png)

add_subdirectory(game_of_life)
