####################

if(WITH_METEOIO)

  find_path(METEOIO_INCLUDE_DIR meteoio/Config.h
    PATHS
    /usr/include
    /usr/local/include
    ${METEOIO_DIR}/include
    )

  find_library(METEOIO_LIBRARY
    NAMES meteoio
    PATHS /usr/lib64 ${METEOIO_DIR}/lib)

endif()  
