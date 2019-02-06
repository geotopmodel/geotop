# Installation instructions

The following instructions have been tested on a Linux system.

## Download the GEOtop v3.0 source code

Clone the git repository:
```
git clone https://github.com/geotopmodel/geotop.git
```

Move to your local geotop repository
```
cd geotop
```

Checkout the v3.0 branch and make sure you are in the right branch:
```
git checkout v3.0
git branch
```


## Install using build system tools
Now you can compile using a build system tool.
Build tools are programs that automate the creation of executable applications
from source code.
Building incorporates compiling, linking and packaging the code into
a usable or executable form.
You can choose between [CMake](https://cmake.org/) (Please note that you need cmake version 3 at least)
and [meson](http://mesonbuild.com/).

## Option 1: install using CMake
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
Now the current directory will have the following files and folders:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/make-build[v3.0*] $ ls -l
totale 84
-rw-rw-r-- 1 elisa elisa 11518 giu 15 14:22 CMakeCache.txt
drwxrwxr-x 5 elisa elisa  4096 giu 15 14:23 CMakeFiles
-rw-rw-r-- 1 elisa elisa  1574 giu 15 14:23 cmake_install.cmake
-rw-rw-r-- 1 elisa elisa   307 giu 15 14:23 CTestTestfile.cmake
-rw-rw-r-- 1 elisa elisa 51860 giu 15 14:23 Makefile
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:22 src
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:23 tests
```
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

## Option 2: Install using Meson
- Create the build directory and go inside it:
```
mkdir meson-build
cd meson-build
```

- Create the build file:
```
meson
```

- Now the current directory will have the following files and folders:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build[v3.0*] $ ls -l
totale 156
-rw-rw-r-- 1 elisa elisa 64560 giu 15 14:39 build.ninja
-rw-rw-r-- 1 elisa elisa 67030 giu 15 14:39 compile_commands.json
drwxrwxr-x 7 elisa elisa  4096 giu 15 14:39 meson
drwxrwxr-x 2 elisa elisa  4096 giu 15 14:39 meson-logs
drwxrwxr-x 2 elisa elisa  4096 giu 15 14:39 meson-private
drwxrwxr-x 4 elisa elisa  4096 giu 15 14:39 src
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:39 subprojects
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:39 tests
```

- Check the default values for the options, opening the file meson.build
in the upper directory or typing:
```
meson configure
```

- If you want to modify some of them, add -Doption=value: for example
    - set the build type to debug type: ``` meson configure -Dbuildtype=debug ```
    - add compiler and linker options (i.e. add ```-pg```): ``` meson configure -Dcpp_args=-pg -Dcpp_link_args=-pg```
    - define multiple compiler options: ```meson configure -Dcpp_args=" -OPTION_1 -OPTION_2"```

- To check if the desired flags were activated you can look at their current values
(```true``` or ```false```) again typing inside the build folder:
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

### Problems
If for some reasons at a certain point after typing ```ninja``` you get a message like:
```
Something went terribly wrong. Please file a bug.
FAILED: build.ninja
```
try to remove the build folder and create it again.
