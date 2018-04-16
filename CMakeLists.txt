cmake_minimum_required(VERSION 2.8)
SET( CMAKE_BUILD_TYPE RELEASE CACHE STRING "A variable which controls the type of build" )
SET( CMAKE_CXX_FLAGS " -std=c++11 " CACHE STRING "")
SET( CMAKE_CXX_FLAGS_DEBUG " -Wextra -Wall -g -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -O0 " CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_RELEASE " -Wextra -Wall -O2 -funroll-loops -fstrict-aliasing " CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_MINSIZEREL " -Wextra -Wall -O2 " CACHE STRING "" )

project(geotop LANGUAGES CXX)

include_directories(src/geotop)
include_directories(src/libraries/ascii)
include_directories(src/libraries/fluidturtle)
include_directories(src/libraries/geomorphology)
include_directories(src/libraries/math)

add_executable(geotop
        src/geotop/blowingsnow.cc
        src/geotop/channels.cc
        src/geotop/clouds.cc
        src/geotop/deallocate.cc
        src/geotop/energy.balance.cc
        src/geotop/geotop.cc
        src/geotop/indices.cc
        src/geotop/input.cc
        src/geotop/meteo.cc
        src/geotop/meteodata.cc
        src/geotop/meteodistr.cc
        src/geotop/output.cc
        src/geotop/parameters.cc
        src/geotop/PBSM.cc
        src/geotop/pedo.funct.cc
        src/geotop/radiation.cc
        src/geotop/recovering.cc
        src/geotop/snow.cc
        src/geotop/tables.cc
        src/geotop/times.cc
        src/geotop/turbulence.cc
        src/geotop/vegetation.cc
        src/geotop/water.balance.cc
        src/libraries/ascii/import_ascii.cc
        src/libraries/ascii/init.cc
        src/libraries/ascii/rw_maps.cc
        src/libraries/ascii/tabs.cc
        src/libraries/ascii/write_ascii.cc
        src/libraries/fluidturtle/alloc.cc
        src/libraries/fluidturtle/datamanipulation.cc
        src/libraries/fluidturtle/error.cc
        src/libraries/fluidturtle/linearalgebra.cc
        src/libraries/fluidturtle/list.cc
        src/libraries/fluidturtle/probability.cc
        src/libraries/fluidturtle/random.cc
        src/libraries/fluidturtle/statistics.cc
        src/libraries/fluidturtle/t_io.cc
        src/libraries/fluidturtle/tensors3D.cc
        src/libraries/fluidturtle/utilities.cc
        src/libraries/fluidturtle/write_dem.cc
        src/libraries/geomorphology/dtm_resolution.cc
        src/libraries/geomorphology/geomorphology.0875.cc
        src/libraries/geomorphology/geomorphology.cc
        src/libraries/geomorphology/networks.cc
        src/libraries/geomorphology/shadows.cc
        src/libraries/math/sparse_matrix.cc
        src/libraries/math/util_math.cc)

