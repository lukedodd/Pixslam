Pixslam
======= 

Pixslam is a toy language for image processing.

Pixslam has two main features. Firstly it is dedicated to image processing, pixslam is designed for elegantly expressing programs which operate pixel wise on images. Secondly pixslam features just in time compilation to x86-64 code, this allows for rapid experimentation with image processing which would otherwise be impossibly slow in an interpreted language.

Right now Pixslam is very much a toy. While more features are planned it is important to note this project was made for my own entertainment rather than to be _useful_.

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

For now this project is documented by examples. The build process runs lots of example Pixslam code. After the build described above you should see an `examples` directory in your build directory. This will be filled with `.psm` Pixslam source files, generated images, and shell or batch files which show how Pixslam was called to generate each output.
