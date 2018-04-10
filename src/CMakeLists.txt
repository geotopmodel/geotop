PROJECT(geotop)

###########################################################
# set variables for internal projects

# ascii
SET( ASCII_LIBRARIES
    "ascii" )
	SET(ASCII_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/libraries/ascii")
	SET(ASCII_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/libraries/ascii")

# gt_utilities
SET( GT_UTILITIES_LIBRARIES
    "gt_utilities" )
	SET( GT_UTILITIES_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/gt_utilities")
	SET( GT_UTILITIES_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/gt_utilities")

# meteoio_plugin
SET( METEOIO_PLUGIN_LIBRARIES 
    "meteoioplugin" )
	SET( METEOIO_PLUGIN_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/meteoio_plugin")
	SET( METEOIO_PLUGIN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/meteoio_plugin")

# geotop
SET( GEOTOP_LIBRARIES
    "geotopAll" )
	SET( GEOTOP_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/src/geotop")
	SET( GEOTOP_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/src/geotop")

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
	${GEOTOP_INCLUDE_PATH}
	${GT_UTILITIES_INCLUDE_PATH}
	${METEOIO_INCLUDE_PATH}
	${METEOIO_PLUGIN_INCLUDE_PATH}
	${Boost_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Energy
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Glacier
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Meteo
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Snow
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Soil
	${CMAKE_CURRENT_SOURCE_DIR}/geotop/Water)

link_directories(
        "/opt/local/lib"
        "/usr/local/lib"
	${ASCII_LIBRARY_PATH}
	${GEOTOP_LIBRARY_PATH}
	${GT_UTILITIES_LIBRARY_PATH}
	${METEOIO_LIBRARY_PATH}
	${METEOIO_PLUGIN_LIBRARY_PATH}
)

SUBDIRS(
	geotop
	gt_utilities
	libraries
	meteoio_plugin
)
