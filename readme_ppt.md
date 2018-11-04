# Using build system tools
Now you can compile using a build system tool.
Build tools are programs that automate the creation of executable applications
from source code.
Building incorporates compiling, linking and packaging the code into
a usable or executable form.
You can choose between [CMake](https://cmake.org/) (Please note that you need cmake version 3 at least)
and [meson](http://mesonbuild.com/).

## CMake
- Create the build directory and go inside it:
```
mkdir cmake-build
cd cmake-build
```
- Check the default values for the options, opening the file CMakeList.txt
in the upper directory or writing:
```
ccmake ..
```
- Press [c] and [e] to configure and edit the options

- Press [t] to toggle the advanced mode; several options will appear.
Select the build type, writing RELEASE or DEBUG after ```CMAKE_BUILD_TYPE```,
and modify the other flags as you prefer, knowing that a flag like:
    - **RELEASE**: will be applied only when compiling in RELEASE mode
    - **DEBUG**: will be applied only when compiling in DEBUG mode.

- Press again [c] and [e] to configure; then press [g] to generate and exit.

- Compile (-j4 allows the usage of 4 processes):
```
make -j4
```

### Testing
- Know which tests are available:
```
 ctest -N
 ```

- Run a single test (i.e. Mazia):
```
ctest -R Mazia
```

- Run a group of tests (i.e. all 1D tests, using 4 processes):
```
ctest -R 1D -j4
```

- Run all tests
```
ctest
```

## Meson
- Create the build directory and go inside it:
```
mkdir meson-build
cd meson-build
```

- Create the build file:
```
meson
```

- Check the default values for the options, opening the file meson.build
in the upper directory or typing:
```
meson configure
```

- Compile:
```
ninja
```

### Testing
- Know which tests are available:
```
 meson test --list
 ```

- Run a single test (i.e. Mazia):
```
meson test --suite geotop:Mazia
```

- Run a group of tests (i.e. all 1D tests):
```
meson test --suite geotop:1D
```

- Run all tests
```
ninja test
```
