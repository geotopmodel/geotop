#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#include "../geotop/constants.h"
#include "../geotop/geotop_common.h"
#include "../geotop/part_snow.h"
#include "../geotop/inputKeywords.h"
#include "../gt_utilities/read_command_line.h"

#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>
#include <iostream>
#include <vector>

void meteoio_init(mio::IOManager& io);
void meteoio_deallocate();

void meteoio_readDEM(GeoMatrix<double>& matrix);
void meteoio_readMap(const std::string& filename, GeoMatrix<double>& matrix);
void meteoio_read2DGrid(GeoMatrix<double>& myGrid, char* _filename);

void meteoio_writeEsriasciiMap(const std::string& filename, GeoMatrix<double>& gm);
void meteoio_writeEsriasciiVector(const std::string& filenam, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV);

void hnw_correction(Par* par, std::vector<mio::MeteoData>& meteo);
void merge_meteo_data(mio::Date& current, std::vector<mio::MeteoData>& meteo);
void meteoio_interpolate(Par* par, double JDbeg, Meteo* met, Water* wat);
void meteoio_interpolate_pointwise(Par* par, double currentdate,	Meteo* met, Water* wat);

void meteoio_interpolate_cloudiness(Par* par, const double& JDbeg, GeoMatrix<double>& tau_cloud_grid, GeoVector<double>& tau_cloud_vec);

bool iswr_present(const std::vector<mio::MeteoData>& vec_meteo, const bool& first_check, AllData *A);

void copyGridToMatrix(mio::Grid2DObject& gridObject, GeoMatrix<double>& myGrid);
void copyGridToMatrix(mio::Grid2DObject& gridObject, GeoMatrix<double>& myGrid, double& thr_min, double& val_min, double& thr_max, double& val_max);
void copyGridToMatrixPointWise(const std::vector<double>& pointValues, GeoMatrix<double>& myGrid);
void copyGridToMatrixPointWise(const std::vector<double>& pointValues, GeoMatrix<double>& myGrid, double& thr_min, double& val_min, double& thr_max, double& val_max);
void changeRHgrid(mio::Grid2DObject& g2d);
void convertToCelsius(mio::Grid2DObject& g2d);
void convertToMBar(mio::Grid2DObject& g2d);
void changeVWgrid(mio::Grid2DObject& g2d, const double& vwMin);
void changeGrid(mio::Grid2DObject& g2d, const double& val);
double tDew(double T, double RH, double P);

#endif
