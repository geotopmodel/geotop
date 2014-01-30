#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#include "../geotop/constants.h"
#include "../geotop/geotop_common.h"
#include "../geotop/part_snow.h"

#include <meteoio/MeteoIO.h>
#include <iostream>
#include <vector>

void meteoio_init(mio::IOManager& io);
void meteoio_readDEM(GeoMatrix<double>& matrix);
void meteoio_readMap(const std::string& filename, GeoMatrix<double>& matrix);
void meteoio_read2DGrid(TInit* pUV, GeoMatrix<double>& myGrid, char* _filename);

void meteoio_writeEsriasciiMap(const std::string& filename, TInit* pUV, GeoMatrix<double>& gm, long pNumber_novalue);
void meteoio_writeEsriasciiVector(const std::string& filenam, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue);

void hnw_correction(Par* par, std::vector<mio::MeteoData>& meteo);
void meteoio_interpolate(Par* par, double JDbeg, Meteo* met, Water* wat);
void meteoio_interpolate_pointwise(Par* par, double currentdate,	Meteo* met, Water* wat);

void meteoio_interpolate_cloudiness(TInit* pUV, Par* par, const double& JDbeg, GeoMatrix<double>& tau_cloud_grid, GeoVector<double>& tau_cloud_vec);

bool iswr_present(const std::vector<mio::MeteoData>& vec_meteo, const bool& first_check, AllData *A);

void copyGridToMatrix(mio::Grid2DObject& gridObject, GeoMatrix<double>& myGrid);
void copyGridToMatrixPointWise(const std::vector<double>& pointValues, GeoMatrix<double>& myGrid);
void changeRHgrid(mio::Grid2DObject& g2d);
void changeTAgrid(mio::Grid2DObject& g2d);
void changePgrid(mio::Grid2DObject& g2d);
void changeVWgrid(mio::Grid2DObject& g2d, double vwMin);
void changeGrid(mio::Grid2DObject& g2d, const double val);
double tDew(double T, double RH, double P);

#endif
