#ifndef __METEOIOPLUGIN_H__
#define __METEOIOPLUGIN_H__

#ifdef __cplusplus
#include <MeteoIO.h> /* C++ header */
#include <iostream>  /* C++ header */
#include <vector>
extern "C" {
#endif

#include "../libraries/ascii/tabs.h"
#include "../geotop/constants.h"
#include "../geotop/struct.geotop.h"
#include "../geotop/meteo.h"

#define PI 3.14159265358979

	void meteoio_init();
	double ***meteoio_readMeteoData(long*** column, METEO_STATIONS *stations, 
							  double novalue, long nrOfVariables, PAR *par, TIMES *times);
	DOUBLEMATRIX *meteoio_readDEM(T_INIT** UV);
	DOUBLEMATRIX *meteoio_read2DGrid(T_INIT* UV, char* _filename);
	void meteoio_interpolate( PAR* par, double JDbeg, METEO* met, WATER* wat);
	void meteoio_interpolate_pointwise(PAR* par, double currentdate,
			METEO* met, WATER* wat);
	void meteoio_interpolate_cloudiness(T_INIT* UV, PAR* par, double JDbeg, DOUBLEMATRIX* tau_cloud_grid, DOUBLEVECTOR* tau_cloud_vec);
#ifdef __cplusplus
	void initializeMetaData(const std::vector<mio::StationData>& vecStation, 
					    const mio::Date& startDate, const double& novalue, PAR *par, METEO_STATIONS *stations);

	void copyGridToMatrix(mio::Grid2DObject& gridObject, DOUBLEMATRIX* myGrid);
	void copyGridToMatrixPointWise(std::vector<double>& pointValues, DOUBLEMATRIX* myGrid);
	void changeRHgrid(mio::Grid2DObject& g2d);
	void changeTAgrid(mio::Grid2DObject& g2d);
	void changePgrid(mio::Grid2DObject& g2d);
	void changeVWgrid(mio::Grid2DObject& g2d, double vwMin);


}
#endif

#endif
