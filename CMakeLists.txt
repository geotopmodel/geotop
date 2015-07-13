cmake_minimum_required(VERSION 2.6)

#default compiler flags
SET( CMAKE_BUILD_TYPE RELEASE CACHE STRING "A variable which controls the type of build" )
SET( CMAKE_CXX_FLAGS_RELEASE "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_C_FLAGS_RELEASE "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_DEBUG "-Wextra -Wall -g3 -O0" CACHE STRING "" )
SET( CMAKE_C_FLAGS_DEBUG "-Wextra -Wall -g3 -O0" CACHE STRING "" )

project(geotop)

# Version of library
SET(VERSION_MAJOR 2 CACHE INTERNAL "")
SET(VERSION_MINOR 0 CACHE INTERNAL "")
SET(VERSION_PATCH 1 CACHE INTERNAL "")

# config.h variables
SET(PACKAGE "geotop")
SET(PACKAGE_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE} ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(PACKAGE_BUGREPORT "stefano.end@gmail.com")
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config_cmake.h.in ${PROJECT_SOURCE_DIR}/config.h)

###########################################################
# user options

SET(METEOIO_PATH "/usr/local" CACHE PATH "location of the MeteoIO library" )
SET(BUILD_STATIC OFF CACHE BOOL "choose between static or dynamic linking")
SET(PRINT_DOUBLE_PRECISION OFF CACHE BOOL "enable printing with 12 significant cyphers")
SET(VERBOSE OFF CACHE BOOL "enable printing of verification data")
SET(STAGED_FOR_REMOVING OFF CACHE BOOL "enable code staged for removing")
SET(ENABLE_INTERNAL_METEODISTR ON CACHE BOOL "Use internal meteo interpolation routine instead of using MeteoIO")
SET(USE_NETCDF OFF CACHE BOOL "enable the use of NETCDF")
SET(METEOIO_OUTPUT OFF CACHE BOOL "enable the use of MeteoIO for output (Experimental)")

###########################################################
# libraries linking type and extension


IF (BUILD_STATIC)
	SET(LINKING_MODE STATIC)
ELSE (BUILD_STATIC)
	SET(LINKING_MODE SHARED)
ENDIF (BUILD_STATIC)

MESSAGE("Linking mode: ${LINKING_MODE}")

#The following pre-processor variables are automatically defined:
# __linux, __GNUC__, __CYGWIN__, _WIN32, __MINGW32__, __APPLE__
IF(UNIX)
	IF(APPLE)
		SET(DLL_EXT "dylib*")
		SET(STAT_EXT "a")
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -Wno-sign-conversion") #for now, keep clang quiet
		SET(CMAKE_SHARED_LINKER_FLAGS "-framework CoreFoundation")

		#set a full rpath for both build & install tree
		SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
		SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
		SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

		# add the automatically determined parts of the RPATH
		# which point to directories outside the build tree to the install RPATH
		SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

		# the RPATH to be used when installing, but only if it's not a system directory
		LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		IF("${isSystemDir}" STREQUAL "-1")
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		ENDIF("${isSystemDir}" STREQUAL "-1")
	ELSE(APPLE)
		SET(DLL_EXT "so*")
		SET(STAT_EXT "a")
	ENDIF(APPLE)
ELSE(UNIX)
	IF(WIN32)
		SET(DLL_EXT "dll")
		SET(STAT_EXT "a")
	ENDIF(WIN32)
ENDIF(UNIX)

###########################################################
# find dependecies

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

###########################################################
# set flags

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

IF (METEOIO_OUTPUT)
    add_definitions(-DMETEOIO_OUTPUT)
ENDIF (METEOIO_OUTPUT)

SUBDIRS(src)

###########################################################
#for the install target
FILE(GLOB hGeotop "src/geotop/*.h")
FILE(GLOB hGeotopEnergy "src/geotop/Energy/*.h")
FILE(GLOB hGeotopGlacier "src/geotop/Glacier/*.h")
FILE(GLOB hGeotopMeteo "src/geotop/Meteo/*.h")
FILE(GLOB hGeotopSnow "src/geotop/Snow/*.h")
FILE(GLOB hGeotopSoil "src/geotop/Soil/*.h")
FILE(GLOB hGeotopWater "src/geotop/Water/*.h")
FILE(GLOB hGt_utilities "src/gt_utilities/*.h")
FILE(GLOB hLibrariesAscii "src/libraries/ascii/*.h")
FILE(GLOB hLibrariesFT "src/libraries/fluidturtle/*.h")
FILE(GLOB hMeteoio_plugin "src/meteoio_plugin/*.h")
FILE(GLOB hNetCDF "src/netCDF/*.h")

INSTALL(FILES ${hGeotop} DESTINATION include/geotop/src/geotop COMPONENT headers)
INSTALL(FILES ${hGeotopEnergy} DESTINATION include/geotop/src/geotop/Energy COMPONENT headers)
INSTALL(FILES ${hGeotopGlacier} DESTINATION include/geotop/src/geotop/Glacier COMPONENT headers)
INSTALL(FILES ${hGeotopMeteo} DESTINATION include/geotop/src/geotop/Meteo COMPONENT headers)
INSTALL(FILES ${hGeotopSnow} DESTINATION include/geotop/src/geotop/Snow COMPONENT headers)
INSTALL(FILES ${hGeotopSoil} DESTINATION include/geotop/src/geotop/Soil COMPONENT headers)
INSTALL(FILES ${hGeotopWater} DESTINATION include/geotop/src/geotop/Water COMPONENT headers)
INSTALL(FILES ${hGt_utilities} DESTINATION include/geotop/src/gt_utilities COMPONENT headers)
INSTALL(FILES ${hLibrariesAscii} DESTINATION include/geotop/src/libraries/ascii COMPONENT headers)
INSTALL(FILES ${hLibrariesFT} DESTINATION include/geotop/src/libraries/fluidturtle COMPONENT headers)
INSTALL(FILES ${hMeteoio_plugin} DESTINATION include/geotop/src/meteoio_plugin COMPONENT headers)
INSTALL(FILES ${hNetCDF} DESTINATION include/geotop/src/netCDF COMPONENT headers)

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove CMakeCache.txt CPackSourceConfig.cmake CPackConfig.cmake CTestTestfile.cmake DartConfiguration.tcl cmake_install.cmake cmake_uninstall.cmake install_manifest.txt Makefile
	COMMAND cmake -E remove "lib/libgeotop*.${DLL_EXT}" "lib/libgeotop*.${STAT_EXT}" "lib/libgeotop*.lib"
	COMMAND cmake -E remove */CMakeCache.txt */CTestTestfile.cmake */cmake_install.cmake
	COMMAND cmake -E remove */*/CMakeCache.txt */*/CTestTestfile.cmake */*/cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	COMMAND cmake -E remove_directory src/CMakeFiles
	COMMAND cmake -E remove_directory src/app/CMakeFiles
	COMMAND cmake -E remove_directory src/app/test/CMakeFiles
	COMMAND cmake -E remove_directory src/app/test/testInputKeywords/CMakeFiles
	COMMAND cmake -E remove_directory src/app/test/unit_test/CMakeFiles
	COMMAND cmake -E remove_directory src/geotop/CMakeFiles
	COMMAND cmake -E remove_directory src/gt_utilities/CMakeFiles
	COMMAND cmake -E remove_directory src/libraries/CMakeFiles
	COMMAND cmake -E remove_directory src/libraries/ascii/CMakeFiles
	COMMAND cmake -E remove_directory src/libraries/fluidturtle/CMakeFiles
	COMMAND cmake -E remove_directory src/meteoio_plugin/CMakeFiles
	COMMAND cmake -E remove_directory src/netCDF/CMakeFiles
)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
