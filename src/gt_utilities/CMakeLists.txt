PROJECT(gt_utilities)

IF(COMMAND CMAKE_POLICY)
        CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

ADD_LIBRARY( gt_utilities ${LINKING_MODE}
			gt_netcdf_utilities.h
			gt_symbols.h
			gt_utilities.h
			ncgt_netcdf2turtle.cc
			ncgt_netcdf2turtle.h
			ncgt_output.cc
			ncgt_output.h
			ncgt_turtle2netcdf.cc
			ncgt_turtle2netcdf.h
			ncgt_utilities.cc
			ncgt_utilities.h
			read_command_line.cc
			read_command_line.h
			read_command_line_netcdf.cc
			read_command_line_netcdf.h
			path_utils.h
			path_utils.c
)

TARGET_LINK_LIBRARIES( gt_utilities
)

SET_TARGET_PROPERTIES(gt_utilities PROPERTIES
  MACOSX_RPATH "${CMAKE_INSTALL_RPATH}"
  PREFIX "${LIBPREFIX}"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib/"
  CLEAN_DIRECT_OUTPUT 1
  VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
  SOVERSION "${VERSION_MAJOR}"
  )
INSTALL(TARGETS gt_utilities
  RUNTIME DESTINATION lib/ COMPONENT libraries
  LIBRARY DESTINATION lib/ COMPONENT libraries
  ARCHIVE DESTINATION lib/ COMPONENT libraries
  )
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

