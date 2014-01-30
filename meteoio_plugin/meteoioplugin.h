#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#include <meteoio/MeteoIO.h> /* C++ header */
#include <iostream>  /* C++ header */
#include <vector>

#include "../geotop/constants.h"
#include "../geotop/geotop_common.h"
#include "../geotop/part_snow.h"

#define PI 3.14159265358979

	void meteoio_init(mio::IOManager& io);

	//	double ***meteoio_readMeteoData(long*** column, METEO_STATIONS *stations,
    //							   long nrOfVariables, PAR *par, TIMES *times);
	double ***meteoio_readMeteoData(long*** column, MeteoStations *stations,
							   long nrOfVariables, Par *par, Times *times);

//	DOUBLEMATRIX *meteoio_readDEM(T_INIT** UV);
	DOUBLEMATRIX *meteoio_readDEM(TInit** UV);
	void meteoio_readDEM(DOUBLEMATRIX*& matrix);
	void meteoio_readDEM(GeoMatrix<double>& matrix);

//	void meteoio_readMap(const std::string& filename, DOUBLEMATRIX*& matrix);
	void meteoio_readMap(const std::string& filename, GeoMatrix<double>& matrix);

//	DOUBLEMATRIX *meteoio_read2DGrid(T_INIT* pUV, char* _filename);
	void meteoio_read2DGrid(TInit* pUV, GeoMatrix<double>& myGrid, char* _filename);

	void meteoio_writeEsriasciiMap(const std::string& filename, TInit* pUV, GeoMatrix<double>& gm, long pNumber_novalue);
	void meteoio_writeEsriasciiVector(const std::string& filenam, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue);

	void hnw_correction(Par* par, std::vector<mio::MeteoData>& meteo);
//	void meteoio_interpolate( PAR* par, double JDbeg, METEO* met, WATER* wat);
	void meteoio_interpolate( Par* par, double JDbeg, Meteo* met, Water* wat);

//	void meteoio_interpolate_pointwise(PAR* par, double currentdate,
//			METEO* met, WATER* wat);
	void meteoio_interpolate_pointwise(Par* par, double currentdate,
			Meteo* met, Water* wat);

//	void meteoio_interpolate_cloudiness(T_INIT* pUV, PAR* par, double JDbeg, DOUBLEMATRIX* tau_cloud_grid, DOUBLEVECTOR* tau_cloud_vec);
	void meteoio_interpolate_cloudiness(TInit* pUV, Par* par, double JDbeg, GeoMatrix<double>& tau_cloud_grid, GeoVector<double>& tau_cloud_vec);

//	void initializeMetaData(const std::vector<mio::StationData>& vecStation,
//					    const mio::Date& startDate, PAR *par, METEO_STATIONS *stations);
	void initializeMetaData(const std::vector<mio::StationData>& vecStation, 
					    const mio::Date& startDate, Par *par, MeteoStations *stations);

	bool iswr_present(const std::vector<mio::MeteoData>& vec_meteo, const bool& first_check, AllData *A);
	void replace_grid_values(const mio::DEMObject& dem, const double& value, mio::Grid2DObject& grid);

	void copyGridToMatrix(mio::Grid2DObject& gridObject, DOUBLEMATRIX* myGrid);
	void copyGridToMatrix(mio::Grid2DObject& gridObject, GeoMatrix<double>& myGrid);
	void copyGridToMatrixPointWise(std::vector<double>& pointValues, DOUBLEMATRIX* myGrid);
	void copyGridToMatrixPointWise(std::vector<double>& pointValues, GeoMatrix<double>& myGrid);
	void changeRHgrid(mio::Grid2DObject& g2d);
	void changeTAgrid(mio::Grid2DObject& g2d);
	void changePgrid(mio::Grid2DObject& g2d);
	void changeVWgrid(mio::Grid2DObject& g2d, double vwMin);
	void changeGrid(mio::Grid2DObject& g2d, const double val);
	void copyInterpMeteoData(double *out, std::vector<mio::MeteoData>& meteoin);
	double tDew(double T, double RH, double P);
	double checkNOvalue(double var);

#endif
