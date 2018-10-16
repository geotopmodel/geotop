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
    -  set the build type to debug type:
```
meson configure -Dbuildtype=debug
```
   - add compiler and linker options (i.e., add ```-pg```)
```
 meson configure -Dcpp_args=-pg -Dcpp_link_args=-pg
```

- Compile:
```
ninja
```
