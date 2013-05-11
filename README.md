Pixslam
======= 

Pixslam is a toy language for image processing.

Pixslam has two main features. 
 * Pixslam is dedicated to image processing, pixslam is designed for elegantly expressing programs which operate pixel wise on images.
 * Pixslam features just in time compilation (JIT) to x86-64 code, this allows for rapid experimentation with image processing which would otherwise be impossibly slow in an interpreted language.

Right now Pixslam is very much a toy. While more features are planned it is important to note this project was made for my own entertainment rather than to be _useful_.

Installation
------------

Standard out of source CMake build. Works on Linux (probably other x86-64 unix OSs too) with recent GCC versions and Windows with Visual Studio 2012. You must build on a x86-64 (otherwise known as amd64) platform!

For Linux builds. Clone the Pixslam repository to somewhere and then do the following:

```bash
mkidr pixslam_build
cd pixslam_build
cmake /path/to/pixlsam/source/folder
make
```

For Windows builds:

 * Clone Pixslam repoistory
 * Fire up cmake gui (install cmake first!) and point it to the Pixslam repo and whatever build directory you want.
 * Click configure, select visual studio 11 win64 project, and then generate.
 * Go to the build directory, open up the pixslam solution file and build all.

Documentation
-------------

For now this project is documented by examples. The build process runs lots of example Pixslam code. After the build described above you should see an `examples` directory in your build directory. This will be filled with `.psm` Pixslam source files, generated images, and shell or batch files which show how Pixslam was run to generate each output.
