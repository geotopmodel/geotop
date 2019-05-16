#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#include "../geotop/constants.h"
#include <meteoio/MeteoIO.h>
#include <iostream>
#include <vector>
#include "matrix.h"

void copyGridToMatrix(mio::Grid2DObject& gridObject, Matrix<double>* matrix); // copy map from MeteoIO to GEOtop

void meteoio_initUV(mio::DEMObject& dem, Matrix<double>* matrix); // copy DEM map from MeteoIO to GEOtop

#endif
