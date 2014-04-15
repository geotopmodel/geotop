PROJECT(fluidturtle)

IF(COMMAND CMAKE_POLICY)
        CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

ADD_LIBRARY( fluidturtle ${LINKING_MODE}
		t_alloc.cc
		t_alloc.h
		t_io.cc
		t_io.h
		tensors3D.cc
		tensors3D.h
		turtle.h
)

TARGET_LINK_LIBRARIES( fluidturtle
		${METEOIO_LIBRARIES}
)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

