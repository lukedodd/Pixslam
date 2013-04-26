
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

ADD_CUSTOM_TARGET(exampeles ALL DEPENDS min_out.png metaballs_out.png)

