cmake_minimum_required(VERSION 2.8)

###########################################################
# general settings

#default compiler flags
SET( CMAKE_BUILD_TYPE RELEASE CACHE STRING "A variable which controls the type of build" )
SET( CMAKE_CXX_FLAGS_DEBUG "-Wextra -Wall -ansi -pedantic -g -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs  -O0" CACHE STRING "" )
SET( CMAKE_C_FLAGS_RELEASE "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_RELEASE "-Wextra -Wall -O2 " CACHE STRING "" )
SET( CMAKE_C_FLAGS_DEBUG "-Wextra -Wall -g3 -O0" CACHE STRING "" )
SET( CMAKE_CXX_FLAGS_MINSIZEREL "-Wextra -Wall -O2" CACHE STRING "" )
SET( CMAKE_C_FLAGS_MINSIZEREL "-Wextra -Wall -O2" CACHE STRING "" )

project(geotop)

# Version of library
SET(VERSION_MAJOR 2 CACHE INTERNAL "")
SET(VERSION_MINOR 1 CACHE INTERNAL "")
SET(VERSION_PATCH 0 CACHE INTERNAL "")

# config.h variables
SET(PACKAGE "geotop")
SET(PACKAGE_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(PACKAGE_STRING "${PACKAGE} ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(PACKAGE_BUGREPORT "stefano.end@gmail.com")
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config_cmake.h.in ${PROJECT_SOURCE_DIR}/config.h)


####################
# define _GNU_SOURCE if C compiler is gnu..

if(CMAKE_COMPILER_IS_GNU_CC)
   SET( _GNU_SOURCE ON)
endif()

###########################################################

IF(APPLE)
    EXEC_PROGRAM(uname ARGS -v OUTPUT_VARIABLE DARWIN_VERSION)
    STRING(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
    MESSAGE("DARWIN_VERSION=${DARWIN_VERSION}")
    IF (DARWIN_VERSION STRLESS 16)
        SET(APPLE_OLD_OSX ON)
    ENDIF(DARWIN_VERSION STRLESS 16)
ENDIF(APPLE)

###########################################################
# user options

SET(METEOIO_PATH "/usr/local/" CACHE PATH "location of the MeteoIO library" )
SET(BUILD_STATIC OFF CACHE BOOL "choose between static or dynamic linking")
SET(PRINT_DOUBLE_PRECISION OFF CACHE BOOL "enable printing with 12 significant cyphers")
SET(VERBOSE OFF CACHE BOOL "enable printing of verification data")
SET(VERY_VERBOSE OFF CACHE BOOL "enable printing of verification data")
SET(ENABLE_INTERNAL_METEODISTR ON CACHE BOOL "Use internal meteo interpolation routine instead of using MeteoIO")
# these below to be checked 
SET(METEOIO_OUTPUT OFF CACHE BOOL "enable the use of MeteoIO for output (Experimental)")

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

IF(BUILD_STATIC)
  FIND_PATH(PROJ_INCLUDE_PATH "proj_api.h")
  FIND_LIBRARY(PROJ_LIBRARIES proj)
  INCLUDE_DIRECTORIES(${PROJ_INCLUDE_PATH})
ENDIF(BUILD_STATIC)

###########################################################
# libraries linking type and extension

IF (BUILD_STATIC)
  SET(LINKING_MODE STATIC)
  SET(BUILD_SHARED_LIBS OFF)
ELSE (BUILD_STATIC)
  SET(LINKING_MODE SHARED)
  SET(BUILD_SHARED_LIBS ON)
ENDIF (BUILD_STATIC)

MESSAGE("Linking mode: ${LINKING_MODE}")

# build code version
INCLUDE("${CMAKE_SOURCE_DIR}/tools/cmake/BuildVersion.cmake")
BuildVersion()

IF(MSVC)
	#SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs"  FORCE)
	SET(WARNINGS "/W4 /D_CRT_SECURE_NO_WARNINGS /EHsc") #Za: strict ansi EHsc: handle c++ exceptions
	#SET(EXTRA_WARNINGS "/Wp64") #/Wall
	SET(OPTIM "/O2 /DNDEBUG /MD /DNOSAFECHECKS")
	SET(ARCH_OPTIM "/arch:SSE2")
	SET(ARCH_SAFE "")
	SET(DEBUG "/Z7 /Od /D__DEBUG /MDd")
	SET(_VERSION "/D_VERSION=\\\"${_versionString}\\\"")
	IF(GUI_EXCEPTIONS)
		SET(MSG_BOX "/DMESG_BOX")
	ENDIF(GUI_EXCEPTIONS)
ELSE(MSVC)
	#we consider that all other compilers support "-" options and silently ignore what they don't know
	SET(WARNINGS "-Wall -Wno-long-long")
	SET(PROFILING "-pg -fprofile-arcs")
	SET(DEEP_WARNINGS "-Wunused-value -Wshadow -Wpointer-arith -Wconversion -Winline -Wdisabled-optimization") #-Wfloat-equal -Wpadded
	SET(EXTRA_WARNINGS "-Wextra -pedantic -Weffc++ ${DEEP_WARNINGS}") #-Weffc++
	IF(CMAKE_COMPILER_IS_GNUCXX AND NOT ANDROID)
		SET(EXTRA_WARNINGS "${EXTRA_WARNINGS} -ansi") #Android does not like ansi...
		IF(WIN32) #for gcc on windows
			SET(CMAKE_SHARED_LINKER_FLAGS "--enable-auto-import")
		ENDIF(WIN32)
	ENDIF(CMAKE_COMPILER_IS_GNUCXX AND NOT ANDROID)
	SET(OPTIM "-g -O3 -DNDEBUG -DNOSAFECHECKS") #-flto for gcc>4.5?
	EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	IF (GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
		SET(ARCH_OPTIM  "-march=native -mtune=native")
	ENDIF()
	SET(ARCH_SAFE  "")
	SET(DEBUG "-g3 -O0 -D__DEBUG")
	SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
ENDIF(MSVC)

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
# set extra flags

IF (PRINT_DOUBLE_PRECISION)
	add_definitions(-DUSE_DOUBLE_PRECISION_OUTPUT)
ENDIF (PRINT_DOUBLE_PRECISION)

IF (VERBOSE)
	add_definitions(-DVERBOSE)
ENDIF (VERBOSE)

IF(ENABLE_INTERNAL_METEODISTR)
	add_definitions(-DUSE_INTERNAL_METEODISTR)
ENDIF(ENABLE_INTERNAL_METEODISTR)

IF (WITH_LOGGER)
    add_definitions(-DWITH_LOGGER)
ENDIF (WITH_LOGGER)

IF (METEOIO_OUTPUT)
    add_definitions(-DMETEOIO_OUTPUT)
ENDIF (METEOIO_OUTPUT)

#Add subdirectory, so it can inherit all variables defined above
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
FILE(GLOB hMeteoio_plugin "src/meteoio_plugin/*.h")

INSTALL(FILES ${hGeotop} DESTINATION include/geotop/src/geotop COMPONENT headers)
INSTALL(FILES ${hGeotopEnergy} DESTINATION include/geotop/src/geotop/Energy COMPONENT headers)
INSTALL(FILES ${hGeotopGlacier} DESTINATION include/geotop/src/geotop/Glacier COMPONENT headers)
INSTALL(FILES ${hGeotopMeteo} DESTINATION include/geotop/src/geotop/Meteo COMPONENT headers)
INSTALL(FILES ${hGeotopSnow} DESTINATION include/geotop/src/geotop/Snow COMPONENT headers)
INSTALL(FILES ${hGeotopSoil} DESTINATION include/geotop/src/geotop/Soil COMPONENT headers)
INSTALL(FILES ${hGeotopWater} DESTINATION include/geotop/src/geotop/Water COMPONENT headers)
INSTALL(FILES ${hGt_utilities} DESTINATION include/geotop/src/gt_utilities COMPONENT headers)
INSTALL(FILES ${hLibrariesAscii} DESTINATION include/geotop/src/libraries/ascii COMPONENT headers)
#INSTALL(FILES ${hLibrariesFT} DESTINATION include/geotop/src/libraries/fluidturtle COMPONENT headers)
INSTALL(FILES ${hMeteoio_plugin} DESTINATION include/geotop/src/meteoio_plugin COMPONENT headers)
#INSTALL(FILES ${hNetCDF} DESTINATION include/geotop/src/netCDF COMPONENT headers)

###########################################################
#for the uninstall target

ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake cmake_uninstall.cmake install_manifest.txt Makefile
	COMMAND cmake -E remove "lib/libgeotop*.${DLL_EXT}" "lib/libgeotop*.${STAT_EXT}" "lib/libgeotop*.lib"
	COMMAND cmake -E remove */CMakeCache.txt */CTestTestfile.cmake */cmake_install.cmake */Makefile
	COMMAND cmake -E remove */*/CMakeCache.txt */*/CTestTestfile.cmake */*/cmake_install.cmake */*/Makefile
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
#	COMMAND cmake -E remove_directory src/netCDF/CMakeFiles
)

###########################################################
#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

      ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

###########################################################
## CPACK section
###########################################################
SET(CPACK_PACKAGE_DESCRIPTION "Hydrological cycle model")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GEOtop is a distributed model of the mass and energy balance of the hydrological cycle, which is applicable to simulations in continuum in small catchments. GEOtop deals with the effects of topography on the interaction between energy balance and hydrological cycle with peculiar solutions.")
# SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_PACKAGE_NAME "GEOtop")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "GEOtop")
SET(CPACK_PACKAGE_CONTACT "Francesco Serafin <f.serafin.1@unitn.it>")
SET(CPACK_PACKAGE_VENDOR "Università degli Studi di Trento")
# SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/icons\\\\icon.ico")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/doc/COPYING.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_SOURCE_DIR}/Welcome.txt")
SET(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.txt")

IF (DEST MATCHES "par")
	SET(PACKAGE_EXTRA "-par")
ENDIF(DEST MATCHES "par")
SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")
SET(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
#SET(CPACK_ARCH "i386")
SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}${PACKAGE_EXTRA}-${VERSION}-${CPACK_ARCH}")

#installation by components for formats that support it
#SET(CPACK_COMPONENTS_ALL_GROUPS_IN_ONE_PACKAGE 1)
IF(APPLE)
	SET(CPACK_COMPONENTS_ALL binary)
ELSE(APPLE)
	SET(CPACK_COMPONENTS_ALL libraries headers documentation examples binary)
ENDIF(APPLE)
SET(CPACK_COMPONENT_LIBRARIES_DISPLAY_NAME "GEOtop")
SET(CPACK_COMPONENT_LIBRARIES_DESCRIPTION "GEOtop model")
SET(CPACK_COMPONENT_LIBRARIES_GROUP Runtime)
SET(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "Header files")
SET(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ header files")
SET(CPACK_COMPONENT_HEADERS_GROUP Development)
SET(CPACK_COMPONENT_HEADERS_DEPENDS libraries)
# IF(INSTALL_DOC)
# 	SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
# 	SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION "Documentation")
# 	SET(CPACK_COMPONENT_DOCUMENTATION_GROUP Docs)
# 	SET(CPACK_COMPONENT_EXAMPLES_DISPLAY_NAME "Examples")
# 	SET(CPACK_COMPONENT_EXAMPLES_DESCRIPTION "C++ Source code of various simple examples of applications with their example data")
# 	SET(CPACK_COMPONENT_EXAMPLES_DEPENDS libraries headers)
# 	SET(CPACK_COMPONENT_EXAMPLES_GROUP Docs)
# 	SET(CPACK_COMPONENT_DOCS_GROUP "Documentation")
# 	SET(CPACK_COMPONENT_GROUP_DOCS_DESCRIPTION "How to use GEOtop")
# ENDIF(INSTALL_DOC)

SET(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")
# SET(CPACK_COMPONENT_GROUP_RUNTIME_DESCRIPTION "What you need to run GEOtop")
SET(CPACK_COMPONENT_DEVELOPMENT_GROUP "Development")
# SET(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION "What you need to develop GEOtop")
SET(CPACK_PACKAGE_EXECUTABLES "geotop;GEOtop")

IF(UNIX)
	IF(APPLE)
		SET(CPACK_GENERATOR "PackageMaker")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/opt/local")
	ELSE(APPLE)
		SET(CPACK_GENERATOR "DEB;RPM;TGZ")
		SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
		SET(CPACK_DEBIAN_PACKAGE_NAME "geotop")
		SET(CPACK_DEBIAN_PACKAGE_SECTION "devel")
		IF(BUILD_STATIC)
		  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libproj0, libstdc++6") #objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		ELSE(BUILD_STATIC)
		  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-filesystem-dev, libboost-system-dev, libboost-iostreams-dev, libboost-regex-dev, libboost-program-options-dev, libboost-test-dev") #objdump -p lib/libmeteoio.so | grep NEEDED | tr -s ' ' | cut -d' ' -f3 | tr '\n' ','
		ENDIF(BUILD_STATIC)
		SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${CPACK_ARCH}) #dpkg --print-architecture
		SET(CPACK_RPM_PACKAGE_NAME "GEOtop")
		SET(CPACK_RPM_PACKAGE_LICENSE "LPGLv3")
		SET(CPACK_RPM_PACKAGE_REQUIRES "libstdc++6") #simply libproj for OpenSuse
		SET(CPACK_RPM_PACKAGE_ARCHITECTURE ${CPACK_ARCH})
		SET(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
	ENDIF(APPLE)
ELSE(UNIX)
	# IF(WIN32)
	# ENDIF(WIN32)
ENDIF(UNIX)


SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${VERSION}" CACHE INTERNAL "tarball basename")
SET(CPACK_SOURCE_GENERATOR TGZ)
# The following components are regex's to match anywhere (unless anchored)
# in absolute path + filename to find files or directories to be excluded
# from source tarball.
SET(CPACK_SOURCE_IGNORE_FILES
	"~$"
	"\\\\.svn/"
	"/CMakeFiles/"
	"/CMakeCache*"
	"/cmake_install\\\\.cmake$"
	"/cmake_uninstall\\\\.cmake$"
	"^_CPack_.*/"
	"/CPack*"
	"\\\\.o$"
	"/m4/"

)

INCLUDE(CPack)

###################### CTEST ######################################
### added by S.Cozzini: maybe it should be set in a separate file #
###################################################################

# generate test program
## test_runner.py variables
SET ( GEOTOP_PROGRAM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/geotop" )
SET ( GEOTOP_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
IF(ENABLE_INTERNAL_METEODISTR)
    SET ( METEOIOPREFIX  off )
ELSE(ENABLE_INTERNAL_METEODISTR)
     SET ( METEOIOPREFIX on )
ENDIF(ENABLE_INTERNAL_METEODISTR)
configure_file (${GEOTOP_TEST_DIR}/test_runner.py.cmake.in ${GEOTOP_TEST_DIR}/test_runner.py)



enable_testing()



#this macro finds every subdirectory of the parameter if there is a geotop.inpts file inside, usage: SUBDIRLIST(SUBDIRS ${DIR})
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      IF(EXISTS ${curdir}/${child}/geotop.inpts)
        LIST(APPEND dirlist ${child})
      ENDIF()
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()


SET(TEST_DIRECTORY ${GEOTOP_TEST_DIR})
SET(1D_DIRECTORY ${TEST_DIRECTORY}/1D)
SET(3D_DIRECTORY ${TEST_DIRECTORY}/3D)
SET(TEST_RUNNER ${TEST_DIRECTORY}/test_runner.py)

SUBDIRLIST(1DTESTS ${1D_DIRECTORY})
SUBDIRLIST(3DTESTS ${3D_DIRECTORY})

foreach(TEST ${1DTESTS})
	add_test(NAME "1D.${TEST}"
         WORKING_DIRECTORY ${1D_DIRECTORY}/${TEST}
         COMMAND $<TARGET_FILE:geotop> geotops.inpts )
	
	#testing with test_runner
	add_test(NAME "1D.test_runner.${TEST}"
         WORKING_DIRECTORY ${1D_DIRECTORY}/${TEST}
	 COMMAND nosetests --verbose ${TEST_RUNNER})
 	
 	set_tests_properties(
		"1D.test_runner.${TEST}" PROPERTIES
		DEPENDS "1D.${TEST}")
ENDFOREACH()

foreach(TEST ${3DTESTS})
	add_test(NAME "3D.${TEST}"
         WORKING_DIRECTORY ${3D_DIRECTORY}/${TEST}
         COMMAND $<TARGET_FILE:geotop> geotops.inpts )

	#testing with test_runner
	add_test(NAME "3D.test_runner.${TEST}"
         WORKING_DIRECTORY ${3D_DIRECTORY}/${TEST}
	 COMMAND nosetests --verbose ${TEST_RUNNER})
 	
 	set_tests_properties(
		"3D.test_runner.${TEST}" PROPERTIES
		DEPENDS "3D.${TEST}")
ENDFOREACH()
