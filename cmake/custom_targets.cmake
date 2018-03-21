# Provide "uninstall" target

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

#for the uninstall target
CONFIGURE_FILE(
	"${PROJECT_SOURCE_DIR}/tools/cmake/cmake_uninstall.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	IMMEDIATE @ONLY)

      ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

###########################################################

#
# Provide "indent" target for indenting all headers and source files
#
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./utils/indent
  COMMENT "Indenting all files in the geotop directories"
  )
