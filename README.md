Pixslam
======= 

Pixslam is a toy language for image processing.

Pixslam has three main features. 

* Pixslam is dedicated to image processing, pixslam is designed for elegantly expressing programs which operate pixel wise on images.
* Pixslam features just in time compilation (JIT) to x86-64 code, this allows for rapid experimentation with image processing which would otherwise be impossibly slow in an interpreted language.
* Pixslam is lightweight and has no complex build dependencies. So it could be integrated into other applications easily.

Right now Pixslam is very much a toy. While more features are planned it is important to note this project was made for my own entertainment rather than to be _useful_.

TODO: Link to blog post for more in depth description.

Building Pixslam
----------------

Standard out of source CMake build. Works on Linux (probably other x86-64 unix OSs too) with recent GCC versions and Windows with Visual Studio 2012. You must build on a x86-64 (otherwise known as amd64) platform!

For Linux builds. Clone the Pixslam repository to somewhere and then do the following:

```bash
mkdir pixslam_build
cd pixslam_build
cmake /path/to/pixlsam/source/folder
make
```

For Windows builds:

 * Clone Pixslam repoistory
 * Fire up cmake gui (install cmake first!) and point it to the Pixslam repo and whatever build directory you want.
 * Click configure, select visual studio 11 Win64 project, and then generate.
 * Go to the build directory, open up the pixslam solution file and build all.

Running Pixslam
--------------- 

The build process described above will have generated a `pixslam` executable for you. The `pixslam` command makes it very easy to run Pixslam functions on specified input images and write the resulting image to a file.

Usage:

```bash
pixslam <code> [input-images] <output>
```

Code can either be supplied directly on the command line (don't forget to surround with double quotes), or as a file path to read in. The number of input images read is dependent on the supplied code. The output argument is optional - it defaults to out.png if not specified. 

Help and simple examples for running `pixslam` are displayed if you run it with no arguments.

A reasonable number of input image formats are supported: jpeg, PNG and BMP. Colour input images can be used but Pixslam will treat them as greyscale. The output image type will be PNG regardless of the extension you give it.

Language Description and Examples
---------------------------------

Pixslam uses a lisp style s-expressions. If you have had any experience with these the language should be quite simple to learn.

Code for Pixslam must be a function definition. This consists of a list containing a list of arguments and then an expression which is evaluated for every pixel of the input images: `((Arg1 Arg2 ...) (Expression)`. All input images must be the same size, and the output image is the same size as the input image.

The expression part of a pixslam function defines it's behaviour. It is basically a calculator with support for indexing images. An example should hopefully make this clear.

```
; Add two images together.
((A B)
     ( * 0.5 ; normalise back to [0,1]
        ( + A B) ; add the current pixel from each input image
     )
)
```

The main thing to note about this example is that it shows one method of indexing images: writing an input name an atom (taking no arguments) evaluates to the current pixel int hat image. Pixslam will evaluate the expression for every pixel in your input images and write the result into your output image.

If we run this example after building pixslam by issuing the following command in the build directory:

```
./pixslam examples/compose.psm example_data/lena.png example_data/duck.png
```

The output will be written to `out.png`. Below you can see the inputs and results.

![Add images example](readme_images/img.jpg "Adding two images: left and middle are two input images, right is the result.")


More Information
-----------------

For now this project is documented by examples. Some of those are shown above. The build process runs lots of example Pixslam code. After you build Pixslam you should see an `examples` directory in your build directory. This will be filled with `.psm` Pixslam source files, generated images, and shell or batch files which show how Pixslam was run to generate each output.

