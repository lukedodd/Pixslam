# To create a gif of the png sequence you can run the following command
# in the directory with the generated images.
# (you need imagemagic installed!)
convert -coalesce -resize 256x256 -filter point -delay 5 -loop 0  *.png game_of_life.gif
