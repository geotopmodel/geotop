cmake_minimum_required(VERSION 3.0)

# support policy changes up to version 3.1.0 to suppress some
# unnecessary warnings
if( "${CMAKE_VERSION}" VERSION_EQUAL "3.1" OR
    "${CMAKE_VERSION}" VERSION_GREATER "3.1" )
  CMAKE_POLICY(VERSION 3.1.0)
endif()


SET( CMAKE_BUILD_TYPE RELEASE CACHE STRING "A variable which controls the type of build" )
SET( CMAKE_CXX_FLAGS " -std=c++11 " CACHE STRING "")
SET( CMAKE_CXX_FLAGS_DEBUG " -Wextra -Wall -g -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -O0 " CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_RELEASE " -Wextra -Wall -O2 -funroll-loops -fstrict-aliasing -DNDEBUG" CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_MINSIZEREL " -Wextra -Wall -O2 " CACHE STRING "" )

subdirs(cmake)


project(geotop LANGUAGES CXX C
  VERSION 3.0)


include_directories(src/geotop)
include_directories(src/libraries/ascii)
include_directories(src/libraries/fluidturtle)
include_directories(src/libraries/geomorphology)
include_directories(src/libraries/math)

# of course we need to inform the compiler where to look for
# the config.h file
include_directories(${CMAKE_BINARY_DIR}/src/geotop)


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
        src/geotop/radiation.cc
        src/geotop/recovering.cc
        src/geotop/snow.cc
        src/geotop/struct.geotop.cc
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
        src/libraries/fluidturtle/t_io.cc
        src/libraries/geomorphology/dtm_resolution.cc
        src/libraries/geomorphology/geomorphology.0875.cc
        src/libraries/geomorphology/geomorphology.cc
        src/libraries/math/sparse_matrix.cc
        src/geotop/vegetation.h
        src/geotop/clouds.h
        src/geotop/water.balance.h
        src/geotop/channels.h
        src/geotop/constants.h
        src/geotop/deallocate.h
        src/geotop/blowingsnow.h
        src/geotop/indices.h
        src/geotop/radiation.h
        src/geotop/input.h
        src/geotop/energy.balance.h
        src/geotop/turbulence.h
        src/geotop/meteo.h
        src/geotop/meteodata.h
        src/geotop/struct.geotop.h
        src/geotop/output.h
        src/geotop/meteodistr.h
        src/geotop/warnings.h
        src/geotop/keywords.h
        src/geotop/pedo.funct.h
        src/geotop/snow.h
        src/geotop/times.h
        src/geotop/tables.h
        src/geotop/PBSM.h
        src/geotop/recovering.h
        src/geotop/parameters.h
        src/libraries/math/util_math.h
        src/libraries/math/sparse_matrix.h
        src/libraries/fluidturtle/t_alloc.h
        src/libraries/fluidturtle/t_datamanipulation.h
        src/libraries/fluidturtle/t_io.h
        src/libraries/geomorphology/dtm_resolution.h
        src/libraries/geomorphology/geomorphology.0875.h
        src/libraries/geomorphology/geomorphology.h
        src/libraries/ascii/write_ascii.h
        src/libraries/ascii/import_ascii.h
        src/libraries/ascii/rw_maps.h
        src/libraries/ascii/init.h
        src/libraries/ascii/tabs.h
        src/libraries/ascii/extensions.h
	src/geotop/logger.cc src/geotop/logger.h
	src/geotop/timer.cc src/geotop/timer.h
	src/geotop/vector.h src/geotop/matrix.h src/geotop/tensor.h
	src/geotop/rowview.h src/geotop/matrixview.h
	src/geotop/geotop_asserts.h
	src/geotop/math.optim.h
	${CMAKE_BINARY_DIR}/src/geotop/version.cc)
enable_testing()
subdirs(tests)


subdirs(cmake/final_message)
