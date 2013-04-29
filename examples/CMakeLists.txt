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
   OUTPUT subtract_out.png
   COMMAND pixslam --logCommand subtract.psm metaballs_fancy_out.png box_5x5_out_2.png subtract_out.png
   DEPENDS pixslam subtract.psm metaballs_fancy_out.png box_5x5_out_2.png 
)


add_custom_target(examples ALL DEPENDS copy_examples copy_example_data compose_out.png erode_3x3_out.png dilate_3x3_out.png erode_5x5_out.png dilate_5x5_out.png min_out.png metaballs_out.png metaballs_binary_out.png metaballs_fancy_out.png box_3x3_out_1.png box_3x3_out_2.png box_5x5_out_1.png box_5x5_out_2.png subtract_out.png)



# Special game of life example.
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/game_of_life)

# Loop to generate lots of iterations of game of life.
add_custom_command(
   OUTPUT game_of_life/game_of_life_0.png
   COMMAND pixslam --logCommand game_of_life.psm ../example_data/glider_gun.png game_of_life/game_of_life_0.png
   DEPENDS pixslam game_of_life.psm
)
set( number 0 )
while( number LESS 100 )
    math( EXPR number "${number} + 1" )
    math( EXPR prev "${number} - 1" )
    add_custom_command(
       OUTPUT game_of_life/game_of_life_${number}.png
       COMMAND pixslam --logCommand game_of_life.psm game_of_life/game_of_life_${prev}.png game_of_life/game_of_life_${number}.png
       DEPENDS pixslam game_of_life.psm game_of_life/game_of_life_${prev}.png
    )
endwhile( number LESS 25 )

add_custom_target(example_game_of_life ALL DEPENDS game_of_life/game_of_life_${prev}.png)
