PROJECT(geotop)

###########################################################
# set variables for internal projects

# ascii
SET( ASCII_LIBRARIES
    "ascii" )
IF(APPLE_OLD_OSX)
	SET(ASCII_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/libraries/ascii")
	SET(ASCII_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/libraries/ascii")
ELSE(APPLE_OLD_OSX)
	SET(ASCII_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/libraries/ascii")
	SET(ASCII_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/libraries/ascii")
ENDIF(APPLE_OLD_OSX)

# fluidturtle
SET( FLUIDTURTLE_LIBRARIES
    "fluidturtle" )
IF(APPLE_OLD_OSX)
	SET( FLUIDTURTLE_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/libraries/fluidturtle")
	SET( FLUIDTURTLE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/libraries/fluidturtle")
ELSE(APPLE_OLD_OSX)
	SET( FLUIDTURTLE_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/libraries/fluidturtle")
	SET( FLUIDTURTLE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/libraries/fluidturtle")
ENDIF(APPLE_OLD_OSX)

# gt_utilities
SET( GT_UTILITIES_LIBRARIES
    "gt_utilities" )
IF(APPLE_OLD_OSX)
	SET( GT_UTILITIES_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/gt_utilities")
	SET( GT_UTILITIES_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/gt_utilities")
ELSE(APPLE_OLD_OSX)
	SET( GT_UTILITIES_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/gt_utilities")
	SET( GT_UTILITIES_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/gt_utilities")
ENDIF(APPLE_OLD_OSX)

# meteoio_plugin
SET( METEOIO_PLUGIN_LIBRARIES 
    "meteoioplugin" )
IF(APPLE_OLD_OSX)
	SET( METEOIO_PLUGIN_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/meteoio_plugin")
	SET( METEOIO_PLUGIN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/meteoio_plugin")
ELSE(APPLE_OLD_OSX)
	SET( METEOIO_PLUGIN_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/meteoio_plugin")
	SET( METEOIO_PLUGIN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/meteoio_plugin")
ENDIF(APPLE_OLD_OSX)

# netCDF
SET( NETCDF_LIBRARIES 
    "netcdfio" )
IF(APPLE_OLD_OSX)
	SET( NETCDF_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/netCDF")
	SET( NETCDF_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/netCDF")
ELSE(APPLE_OLD_OSX)
	SET( NETCDF_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/netCDF")
	SET( NETCDF_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/netCDF")
ENDIF(APPLE_OLD_OSX)

# geotop
SET( GEOTOP_LIBRARIES
    "geotopAll" )
IF(APPLE_OLD_OSX)
	SET( GEOTOP_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/geotop")
	SET( GEOTOP_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/geotop")
ELSE(APPLE_OLD_OSX)
	SET( GEOTOP_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/geotop")
	SET( GEOTOP_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/geotop")
ENDIF(APPLE_OLD_OSX)

###########################################################
# set variables for internal libraries

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

###########################################################
# directories

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	${ASCII_INCLUDE_PATH}
	${FLUIDTURTLE_INCLUDE_PATH}
	${GEOTOP_INCLUDE_PATH}
	${GT_UTILITIES_INCLUDE_PATH}
	${METEOIO_INCLUDE_PATH}
	${METEOIO_PLUGIN_INCLUDE_PATH}
	${NETCDF_INCLUDE_PATH}
	${Boost_INCLUDE_DIRS}
) 

link_directories(
        "/opt/local/lib"
        "/usr/local/lib"
	${ASCII_LIBRARY_PATH}
	${FLUIDTURTLE_LIBRARY_PATH}
	${GEOTOP_LIBRARY_PATH}
	${GT_UTILITIES_LIBRARY_PATH}
	${METEOIO_LIBRARY_PATH}
	${METEOIO_PLUGIN_LIBRARY_PATH}
	${NETCDF_LIBRARY_PATH}
)

SUBDIRS(
	geotop
	gt_utilities
	libraries
	meteoio_plugin
	netCDF
)
