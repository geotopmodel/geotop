#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#include "../geotop/constants.h"
#include <meteoio/MeteoIO.h>
#include <iostream>
#include <vector>
#include "matrix.h"

//void meteoio_init(mio::IOManager& iomanager); // (1)

void copyGridToMatrix(mio::Grid2DObject& gridObject, Matrix<double>* matrix); // copy map from MeteoIO to GEOtop

void meteoio_readDEM(mio::DEMObject& dem, Matrix<double>* matrix); // copy DEM map from MeteoIO to GEOtop

//void meteoio_readMap(const std::string &filename, Matrix<double>* matrix); // (4)


#endif
