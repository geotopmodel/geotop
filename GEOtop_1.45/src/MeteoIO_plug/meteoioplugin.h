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


	#define PI 3.14159265358979

	double ***meteoio_readMeteoData(long*** column, METEO_STATIONS *stations, 
							  double novalue, long nrOfVariables, PAR *par, TIMES *times);
	DOUBLEMATRIX *meteoio_readDEM(T_INIT** UV);
	DOUBLEMATRIX *meteoio_read2DGrid(T_INIT* UV, char* _filename);
	DOUBLEMATRIX * doubletens_to_doublemat(DOUBLETENSOR * input,  DOUBLEMATRIX * output, long indx);
	DOUBLETENSOR * doublemat_to_doubletens(DOUBLEMATRIX * input,DOUBLETENSOR * output, long indx);
	void meteoio_interpolate(T_INIT* UV, PAR* par, double JDbeg, METEO* met, WATER* wat, long k);

#ifdef __cplusplus
	void initializeMetaData(const std::vector<mio::StationData>& vecStation, 
					    const mio::Date& startDate, const double& novalue, PAR *par, METEO_STATIONS *stations);

	void copyGridToMatrix(const double& novalue, const mio::Grid2DObject& gridObject, DOUBLEMATRIX* myGrid);
	void copyGridToTensor(const double& novalue, const mio::Grid2DObject& gridObject,DOUBLETENSOR* myGrid, long indx);
	void copyGridToTensorPointWise(const double& novalue, const mio::Grid2DObject& gridObject, DOUBLETENSOR* myGrid, long indx,DOUBLEMATRIX *coordinates);
	void copyGridToMatrixPointWise(const double& novalue, const mio::Grid2DObject& gridObject, DOUBLEMATRIX* myGrid, DOUBLEMATRIX* coordPoints);
	void changeRHgrid(mio::Grid2DObject& g2d);
	void changeTAgrid(mio::Grid2DObject& g2d);
	void changePgrid(mio::Grid2DObject& g2d);


}
#endif

#endif
