####################

if(WITH_METEOIO)
  IF (METEOIO_DIR)
    set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -L${METEOIO_DIR}/lib -I${METEOIO_DIR}/include -lmeteoio" CACHE STRING "" FORCE)
  else()
    set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -lmeteoio" CACHE STRING "" FORCE)  
  endif()
endif()  
