# Special game of life example.

# This is all a little convoluted in camke.
# All we want to do run pixsmal with game_of_life.psm with the glider_gun.png image
# and then run game_of_life.psm on the output from that, and again on the output from that etc...
# 
# Supreme surprise: this ends up simulating Conway's game of life. After a successful build
# you should see a number game_of_life_<number>.png images. Flick though them and you should see a
# glider gun making... gliders!

# Loop to generate lots of iterations of game of life.
copy_files("copy_game_of_life_example" ${CMAKE_CURRENT_SOURCE_DIR}/*.psm ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(
   OUTPUT game_of_life_1001.png
   COMMAND pixslam --logCommand game_of_life.psm ../../example_data/glider_gun.png game_of_life_1001.png
   DEPENDS pixslam game_of_life.psm
)

# Start with a high number so lexicographic ordering does not break.
# (Let's us play in order easily in lots of image viewers to "animate")
set( number 1001 )
while( number LESS 1200 )
    math( EXPR number "${number} + 1" )
    math( EXPR prev "${number} - 1" )
    add_custom_command(
       OUTPUT game_of_life_${number}.png
       COMMAND pixslam --logCommand game_of_life.psm game_of_life_${prev}.png game_of_life_${number}.png
       DEPENDS pixslam game_of_life.psm game_of_life_${prev}.png
    )
endwhile( number LESS 1200 )

add_custom_target(example_game_of_life ALL DEPENDS game_of_life_${prev}.png)
