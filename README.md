Pixslam
======= 

Pixslam is a toy language for image processing.

Pixslam has two main features. Firstly it is dedicated to image processing, pixslam is for writing programs which operate pixel wise on images. Secondly it features just in time compilation to x86-64 code, this allows experimentation with image processing code which can be impossibly slow in a normal dynamic language.

Right now Pixslam is very much a toy. While more features are planned it is important to note this project was made for my my own entertainment and not to be _useful_.

Installation
------------

Standard out of source CMake build. Works on Linux with recent GCC and on Windows with visual studio 2012. You must build on a x86-64 (otherwise known as amd64) platform!

```bash
mkidr pixslam_build
cd pixslam_build
cmake /path/to/pixlsam/source/folder
make
```

Documentation
-------------

Fow now this project is documented by examples. After the build described above you should see an `examples` foldr in your build directory. This will be filled with `.psm` pixslam source files, images, and shell or batch files which show how pixslam was called to generate each output.
