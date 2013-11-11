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
)

TARGET_LINK_LIBRARIES( gt_utilities
)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

