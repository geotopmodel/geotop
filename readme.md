# GEOtop 3.0

## Introduction

GEOtop 3.0 version starts from 2.0 version (branch ```se27xx```) already validated and published
in the Endrizzi et al. 2014 paper.
It performs exactly as the previous version but it has some improvements in terms of:
- usage of object-oriented approach
- development of new data structures
- ease of compiling and running
- modularity and flexibility
- increase in testing coverage.

However, it still lacks of the integration with the MeteoIO library and other features implemented in the current 2.1 version (branch ```master```).
In the next months we plan to move toward a stable 3.0 version, together with a publication.

A more detailed description of this new version can be found in the MHPC thesis (https://www.mhpc.it/)
of Elisa Bortoli at the following link: https://iris.sissa.it/handle/20.500.11767/86154#.XEXzOsZ7l8w

GEOtop 3.0 can be compiled and run following the listed instructions,
tested on a Linux system.

## Getting the source code

- Clone the git repository:
```
git clone https://github.com/geotopmodel/geotop.git
```

- Move to your local geotop repository:
```
cd geotop
```

- Go into the branch v3.0 and make sure you are in:
```
git checkout v3.0
git branch
```

- to update your local repository to the newest commit, execute 
```
git pull origin v3.0
```

## Compiling
Now you can compile using a build system tool.
Build tools are programs that automate the creation of executable applications
from source code.
Building incorporates compiling, linking and packaging the code into
a usable or executable form.
You can choose between [CMake](https://cmake.org/) (Please note that you need cmake version 3 at least)
and [Meson](http://mesonbuild.com/).

### Option 1: using CMake
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
You can modify the values of the flags going to the correspondent line,
pressing "Enter" key and then editing; to save what you have just written
press again "Enter".

- For example you can choose the build type, writing RELEASE (default option) or
DEBUG after ```CMAKE_BUILD_TYPE```,
and modify the other flags as you prefer, knowing that a flag like:
    - *_RELEASE: will be applied only when compiling in RELEASE mode
    - *_DEBUG: will be applied only when compiling in DEBUG mode.

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

### Option 2: using Meson
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
    - to set the build type to debug writes:
    ```
    meson configure -Dbuildtype=debug
    ```
    - to add compiler and linker options (i.e. add ```-pg```) write:
    ```
    meson configure -Dcpp_args=-pg -Dcpp_link_args=-pg
    ```
    - define multiple compiler options:
    ```
    meson configure -Dcpp_args=" -OPTION_1 -OPTION_2"
    ```

- To check if the desired flags were activated, you can look at their current values
(```true``` or ```false```) again typing inside the build folder:
```
meson configure
```

- Compile:
```
ninja
```

## Running the test cases
Now you can run the proposed test cases.

### Option 1: using CMake
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



### Option 2: using Meson
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

## Getting help

### Output error message
If for some reasons at a certain point after typing ```ninja``` you get a message like:
```
Something went terribly wrong. Please file a bug.
FAILED: build.ninja
```
remove the build folder and create it again.

### Looking into the documentation
An interactive documentation can be built with Doxygen
by typing in the root directory:
```
doxygen Doxyfile
```
A new folder **doxygen_generated_doc** will be created, containing two subfolders:
*html* and *latex*, with graphs of the single functions.

If you want to navigate files and functions, go inside
*html* folder and type:
```
firefox index.html
```


### Reporting an issue
To report a problem you can open an issue on GitHub (https://github.com/geotopmodel/geotop/issues) listing all the following infos:
- short description of what happened
- operating system (OS)
- architecture
- compiler
- configuration (type from the build folder ```ccmake ..```
  or ```meson configure```, depending on the build system tool you are using,
  and put the info in a file).
