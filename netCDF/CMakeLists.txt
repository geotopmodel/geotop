PROJECT(netcdfio)

IF(COMMAND CMAKE_POLICY)
        CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

ADD_LIBRARY( netcdfio ${LINKING_MODE}
		gt_symbols.h
		netcdfIO.cc
		netcdfIO.h
		read_command_line.cc
		read_command_line.h
		read_command_line_netcdf.cc
		read_command_line_netcdf.h
)

TARGET_LINK_LIBRARIES( netcdfio
)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

