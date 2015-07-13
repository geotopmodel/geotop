PROJECT(geotop)

# ascii
SET( ASCII_LIBRARIES
    "ascii" )
SET(ASCII_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/libraries/ascii")
SET(ASCII_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/libraries/ascii")

# fluidturtle
SET( FLUIDTURTLE_LIBRARIES
    "fluidturtle" )
SET( FLUIDTURTLE_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/libraries/fluidturtle")
SET( FLUIDTURTLE_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/libraries/fluidturtle")

# gt_utilities
SET( GT_UTILITIES_LIBRARIES
    "gt_utilities" )
SET( GT_UTILITIES_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/gt_utilities")
SET( GT_UTILITIES_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/gt_utilities")

# inputKeywords
SET( INPUT_KEYWORDS_LIBRARIES
    "inputKeywords" )
SET( INPUT_KEYWORDS_LIBRARY_PATH "${PROJECT_LIST_DIR}/src/geotop")
SET( INPUT_KEYWORDS_INCLUDE_PATH "${PROJECT_LIST_DIR}/src/geotop")

# meteoio_plugin
SET( METEOIO_PLUGIN_LIBRARIES 
    "meteoioplugin" )
SET( METEOIO_PLUGIN_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/meteoio_plugin")
SET( METEOIO_PLUGIN_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/meteoio_plugin")

# netCDF
SET( NETCDF_LIBRARIES 
    "netcdfio" )
SET( NETCDF_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/netCDF")
SET( NETCDF_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/netCDF")

# geotop
SET( GEOTOP_LIBRARIES
    "geotopAll" )
SET( GEOTOP_LIBRARY_PATH "${PROJECT_SOURCE_DIR}/src/geotop")
SET( GEOTOP_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/src/geotop")

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

SUBDIRS(
	apps
	geotop
	gt_utilities
	libraries
	meteoio_plugin
	netCDF
)
