cmake_minimum_required(VERSION 2.6)

#default compiler flags
SET( CMAKE_BUILD_TYPE RELEASE CACHE STRING "A variable which controls the type of build" )
SET( CMAKE_CXX_FLAGS_RELEASE "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_C_FLAGS_RELEASE "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_DEBUG "-Wextra -Wall -g3 -O0" CACHE STRING "" )
SET( CMAKE_C_FLAGS_DEBUG "-Wextra -Wall -g3 -O0" CACHE STRING "" )

project("geotop")

# config.h variables
SET(PACKAGE "geotop")
SET(PACKAGE_VERSION "1.2.10")
SET(PACKAGE_STRING "${PACKAGE} 1.2.10")
SET(PACKAGE_BUGREPORT "stefano.end@gmail.com")

SET(METEOIO_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../geotopio/trunk/fakeroot/usr" CACHE PATH "location of the MeteoIO library" )

SET(BUILD_STATIC OFF CACHE BOOL "choose between static or dynamic linking")

SET(PRINT_DOUBLE_PRECISION OFF CACHE BOOL "enable printing with 12 significant cyphers")

SET(VERBOSE OFF CACHE BOOL "enable printing of verification data")

SET(STAGED_FOR_REMOVING OFF CACHE BOOL "enable code staged for removing")

SET(ENABLE_INTERNAL_METEODISTR ON CACHE BOOL "Use internal meteo interpolation, routine instead of using MeteoIO")

SET(USE_NETCDF OFF CACHE BOOL "enable the use of NETCDF")

SET(WITH_LOGGER OFF CACHE BOOL "enable the use of GEOtop Logger (Experimental)")

find_package(PkgConfig REQUIRED)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

IF (BUILD_STATIC)
	SET( LINKING_MODE STATIC)
ELSE (BUILD_STATIC)
	SET( LINKING_MODE SHARED)
ENDIF (BUILD_STATIC)

MESSAGE("Linking mode: ${LINKING_MODE}")

FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(Boost COMPONENTS regex REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
FIND_PACKAGE(Boost COMPONENTS filesystem REQUIRED)
FIND_PACKAGE(Boost COMPONENTS system REQUIRED)
FIND_PACKAGE(Boost COMPONENTS iostreams REQUIRED)
FIND_PACKAGE(Boost COMPONENTS unit_test_framework REQUIRED)

IF(NOT EXISTS "${METEOIO_PATH}/include/meteoio/MeteoIO.h")
	MESSAGE(FATAL_ERROR "METEOIO LIBRARY NOT FOUND")
ELSE(NOT EXISTS "${METEOIO_PATH}/include/meteoio/MeteoIO.h")
	SET(METEOIO_LIBRARY_PATH "${METEOIO_PATH}/lib")
	SET(METEOIO_INCLUDE_PATH "${METEOIO_PATH}/include")
	SET(METEOIO_LIBRARIES "meteoio")
ENDIF(NOT EXISTS "${METEOIO_PATH}/include/meteoio/MeteoIO.h")

IF (PRINT_DOUBLE_PRECISION)
	add_definitions(-DUSE_DOUBLE_PRECISION_OUTPUT)
ENDIF (PRINT_DOUBLE_PRECISION)

IF (VERBOSE)
	add_definitions(-DVERBOSE)
ENDIF (VERBOSE)

IF (STAGED_FOR_REMOVING)
	add_definitions(-DSTAGED_FOR_REMOVING)
ENDIF (STAGED_FOR_REMOVING)

IF(ENABLE_INTERNAL_METEODISTR)
	add_definitions(-DUSE_INTERNAL_METEODISTR)
ENDIF(ENABLE_INTERNAL_METEODISTR)

IF (USE_NETCDF)
	add_definitions(-DUSE_NETCDF)
ENDIF (USE_NETCDF)

IF (WITH_LOGGER)
    add_definitions(-DWITH_LOGGER)
ENDIF (WITH_LOGGER)

#add_definitions(
#    "-std=c++11"
#)

# ascii
SET( ASCII_LIBRARIES
    "ascii" )
SET(ASCII_LIBRARY_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libraries/ascii")
SET(ASCII_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libraries/ascii")

# fluidturtle
SET( FLUIDTURTLE_LIBRARIES
    "fluidturtle" )
SET( FLUIDTURTLE_LIBRARY_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libraries/fluidturtle")
SET( FLUIDTURTLE_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libraries/fluidturtle")

# wrf
SET( WRF_LIBRARIES
    "wrfIO" )
SET( WRF_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/wrf)
SET( WRF_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/wrf)

# gt_utilities
SET( GT_UTILITIES_LIBRARIES
    "gt_utilities" )
SET( GT_UTILITIES_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/gt_utilities)
SET( GT_UTILITIES_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/gt_utilities)

# inputKeywords
SET( INPUT_KEYWORDS_LIBRARIES
    "inputKeywords" )
SET( INPUT_KEYWORDS_LIBRARY_PATH ${CMAKE_CURRENT_LIST_DIR}/geotop)
SET( INPUT_KEYWORDS_INCLUDE_PATH ${CMAKE_CURRENT_LIST_DIR}/geotop)

# meteoio_plugin
SET( METEOIO_PLUGIN_LIBRARIES 
    "meteoioplugin" )
SET( METEOIO_PLUGIN_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/meteoio_plugin)
SET( METEOIO_PLUGIN_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/meteoio_plugin)

# netCDF
SET( NETCDF_LIBRARIES 
    "netcdfio" )
SET( NETCDF_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/netCDF)
SET( NETCDF_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/netCDF)

# geotop
SET( GEOTOP_LIBRARIES
    "geotopAll" )
SET( GEOTOP_LIBRARY_PATH ${CMAKE_CURRENT_SOURCE_DIR}/geotop)
SET( GEOTOP_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/geotop)

# inputKeywords
SET( INPUT_KEYWORDS_LIBRARIES
    "inputKeywords" )

# geotopOutPostProcessing
SET( GEOTOP_OUT_POSTPROC_LIBRARIES
    "geotopOutPostProcessing" )

# geotopOutPostProcessingNew
SET( GEOTOP_OUT_POSTPROC_LIBRARIES_NEW
    "geotopOutPostProcessingNew" )

# geotopMeteo
SET( GEOTOP_METEO_LIBRARIES
    "geotopMeteo" )

# geotopWaterBalance
SET( GEOTOP_WATER_BALANCE_LIBRARIES
    "geotopWaterBalance" )

# geotopPartSnow
SET( GEOTOP_PART_SNOW_LIBRARIES
    "geotopPartSnow" )

# geotopEnergyBalance
SET( GEOTOP_ENERGY_BALANCE_LIBRARIES
    "geotopEnergyBalance" )

# geotopCommon
SET( GEOTOP_COMMON_LIBRARIES
    "geotopCommon" )

# logger
SET( GEOTOP_LOGGER
    "logger" )

# geotopOutputFile
SET( GEOTOP_OUTPUT_FILE
    "geotopOutputFile" )

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${INPUT_KEYWORDS_INCLUDE_PATH}
	${ASCII_INCLUDE_PATH}
	${FLUIDTURTLE_INCLUDE_PATH}
	${GEOTOP_INCLUDE_PATH}
	${GT_UTILITIES_INCLUDE_PATH}
	${METEOIO_INCLUDE_PATH}
	${METEOIO_PLUGIN_INCLUDE_PATH}
	${NETCDF_INCLUDE_PATH}
	${WRF_INCLUDE_PATH}
	${Boost_INCLUDE_DIRS}
) 

link_directories(
        "/opt/local/lib"
        "/usr/local/lib"
	${ASCII_LIBRARY_PATH}
	${INPUT_KEYWORDS_LIBRARY_PATH}
	${FLUIDTURTLE_LIBRARY_PATH}
	${GEOTOP_LIBRARY_PATH}
	${GT_UTILITIES_LIBRARY_PATH}
	${METEOIO_LIBRARY_PATH}
	${METEOIO_PLUGIN_LIBRARY_PATH}
	${NETCDF_LIBRARY_PATH}
	${WRF_LIBRARY_PATH}
)

# generate test program
## test_runner.py variables
SET ( GEOTOP_PROGRAM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/geotop/geotop" )
SET ( GEOTOP_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sample_run")
IF(ENABLE_INTERNAL_METEODISTR)
    SET ( METEODISTR_RESULTS_DIR results.meteodistr )
ELSE(ENABLE_INTERNAL_METEODISTR)
    SET ( METEODISTR_RESULTS_DIR results )
ENDIF(ENABLE_INTERNAL_METEODISTR)
configure_file (${GEOTOP_TEST_DIR}/test_runner.py.cmake.in ${GEOTOP_TEST_DIR}/test_runner.py)
configure_file (${GEOTOP_TEST_DIR}/simple_runner.py.cmake.in ${GEOTOP_TEST_DIR}/simple_runner.py)

# generate config include file
## config.h variables
SET(PACKAGE "geotop")
SET(PACKAGE_VERSION "1.2.10")
SET(PACKAGE_STRING "${PACKAGE} 1.2.10")
SET(PACKAGE_BUGREPORT "stefano.end@gmail.com")
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

SUBDIRS(
	apps
	geotop
	gt_utilities
	libraries
	meteoio_plugin
	netCDF
	wrf
)

