Installation Instructions
=========================

:author: Stefano Cozzini (eXact-lab srl) /Emanuele Cordano 
:date: last revision January 2017


To install GEOtop you need the following libraries and tools:

- Git
- CMake 2.8 or later (optional but recommended: ccmake ncurses ui)
- Boost 1.49 or later (filesystem, system, iostreams, regex,program_options, unit_test_framework and spirit_classic)
- PROJ.4 libraries 4.7.0
- MeteoIO 2.5.1

GEOtop is known to run under the following Operating Systems:

- Mac OS X 10.8 or later
- CentOS 6.5 or later
- Debian 7
- Ubuntu 14.04 or later

CentOS 6.5 and superior
========================= 

1 Install the tools via yum::

    $ su
    # yum install cmake git

2 Download the Boost libraries from www.boost.org

3 Follow the Getting Started Guide on the Boost website to install the libraries

4 Download the PROJ.4 libraries from  http://download.osgeo.org/proj/proj-4.9.3.tar.gz

5 Unpack the tar.gz and run the usual sequence::

    $ tar -xvzf proj-4.9.3.tar.gz
    $ cd proj-4.9.3
    $ ./configure
    $ make
    $ su
    # make install

6 Download the MeteoIO libraries from http://models.slf.ch/p/meteoio/downloads/get/MeteoIO-2.5.1-src.tar.gz

7 Unpack the tar.gz, and run ccmake::
    $ tar -xvzf MeteoIO-2.5.1-src.tar.gz
    $ cd meteoio
    $ ccmake .

8 IMPORTANT: Enable PROJ.4 support in MeteoIO and then configure and
   generate the makefiles.

9 Compile and install MeteoIO::
    $ make
    $ su
    # make install

Now you have installed all the dependecies for GEOtop.

Debian 7
========

1 Install the tools via apt-get::

    $ su
    # apt-get install cmake-curses-gui git

2 Install the Boost libraries via apt-get::

    # apt-get install libboost-all-dev

3 Install the PROJ.4 librareis via apt-get

    # apt-get install libproj-dev

7 Compile and install MeteoIO, see above

Now you have installed all the dependecies for GEOtop.

Compiling GEOtop
================

1 Get the sources from Github::

    $ git clone https://github.com/geotopmodel/geotop.git

2 Run ccmake::

    $ cd geotop
    $ ccmake .

3 Configure and generate the makefiles

4 Run make
5 Run make test to check that all tests work fine 

CMake options
=============

Here the list of option you may specify when you configure Cmake::

     Option                     Default
 BUILD_STATIC                     OFF
 CMAKE_BUILD_TYPE                 RELEASE
 CMAKE_INSTALL_PREFIX             /usr/local
 ENABLE_INTERNAL_METEODISTR       ON
 METEOIO_OUTPUT                   OFF
 METEOIO_PATH                     /usr/local
 VERBOSE                          OFF
 VERY_VERBOSE                     OFF


BUILD_STATIC: If set to ON will build GEOtop as a single static binary.

CMAKE_BUILD_TYPE: Can be RELEASE or DEBUG. The latter sets up GEOtop
logging to high verbosity and enables GDB support.

CMAKE_INSTALL_PREFIX: Sets the "make install" target directory.

ENABLE_INTERNAL_METEODISTR: If set to ON GEOtop will NOT use MeteoIO for
data interpolation.

METEOIO_OUTPUT: If set to ON GEOtop will use MeteoIO to output ARC
files.

METEOIO_PATH: Path to MeteoIO libraries

VERBOSE: increases logging verbosity.

VERY_VERBOSE: increases ieven more logging verbosity.

