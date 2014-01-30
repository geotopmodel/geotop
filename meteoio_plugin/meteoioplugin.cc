#include "meteoioplugin.h"
#include "../geotop/times.h"
#include <sstream>
#include <stdlib.h>
#include "geotop_common.h"
#include "inputKeywords.h"

using namespace std;
using namespace mio;

IOManager* io;
DEMObject dem;

void meteoio_init(mio::IOManager& iomanager) {

	/*
	 *  This function performs the initialization by calling the respective methods within MeteoIO
	 *  1) Reading the DEM (necessary for several spatial interpolations algorithms)
	 */

	io = &iomanager; //pointer to the iomanager instantiated in geotop.cc

	try {
		io->readDEM(dem);
	} catch (exception& e) {
		cerr << "[E] MeteoIO: " << e.what() << endl;
	}
}

void meteoio_readDEM(GeoMatrix<double>& matrix) {
	//copy DEM to topo struct
	matrix.resize(dem.nrows+1, dem.ncols+1);

	geotop::common::Variables::UV->V.resize(2+1);
	geotop::common::Variables::UV->U.resize(4+1);
	geotop::common::Variables::UV->V[1] = -1.0;
	geotop::common::Variables::UV->V[2] = geotop::input::gDoubleNoValue; //GEOtop nodata -9999.0

	geotop::common::Variables::UV->U[1] = dem.cellsize;
	geotop::common::Variables::UV->U[2] = dem.cellsize;
	geotop::common::Variables::UV->U[3] = dem.llcorner.getNorthing();
	geotop::common::Variables::UV->U[4] = dem.llcorner.getEasting();

	copyGridToMatrix(dem, matrix);
}

void meteoio_readMap(const std::string& filename, GeoMatrix<double>& matrix) {
	Grid2DObject temp;
	io->read2DGrid(temp, filename + ".asc");

	matrix.resize(temp.nrows+1, temp.ncols+1);
	copyGridToMatrix(temp, matrix);
}

//DOUBLEMATRIX *meteoio_read2DGrid(T_INIT* pUV, char* _filename) {
void meteoio_read2DGrid(TInit* pUV, GeoMatrix<double>& myGrid, char* _filename) {
//	DOUBLEMATRIX *myGrid = NULL;

	try {

		Grid2DObject gridObject;
		Config cfg("io.ini");
		IOHandler iohandler(cfg);

		std::string gridsrc = "", filename = "";
		cfg.getValue("GRID2D", "INPUT", gridsrc);

		if (gridsrc == "GRASS")
			filename = string(_filename) + ".grassasc";
		else if (gridsrc == "ARC")
			filename = string(_filename) + ".asc";

		iohandler.read2DGrid(gridObject, filename);

	//	if (UV->U->co[1] != gridObject.cellsize)
		if (geotop::common::Variables::UV->U[1] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
	//	else if (UV->U->co[2] != gridObject.cellsize)
		else if (geotop::common::Variables::UV->U[2] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
	//	else if (UV->U->co[3] != gridObject.llcorner.getNorthing())
		else if (geotop::common::Variables::UV->U[3] != gridObject.llcorner.getNorthing())
			throw IOException("Inconsistencies between 2D Grids read", AT);
	//	else if (UV->U->co[4] != gridObject.llcorner.getEasting())
		else if (geotop::common::Variables::UV->U[4] != gridObject.llcorner.getEasting())
			throw IOException("Inconsistencies between 2D Grids read", AT);

		for (unsigned int ii = 0; ii < gridObject.nrows; ii++) {
			for (unsigned int jj = 0; jj < gridObject.ncols; jj++) {
				if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii)
						== IOUtils::nodata) {
				//	myGrid->co[ii + 1][jj + 1] = UV->V->co[2];
					myGrid[ii + 1][jj + 1] = geotop::common::Variables::UV->V[2];
				} else {
				//	myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj,
				//			gridObject.nrows - 1 - ii);
					myGrid[ii + 1][jj + 1] = gridObject.grid2D(jj,
											gridObject.nrows - 1 - ii);
				}
			}
		}
		std::cout << "Read 2D Grid from file '" << filename << "' : "
				<< gridObject.nrows << "(rows) " << gridObject.ncols << "(cols)"
				<< std::endl;
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}
//	return myGrid;
}


void meteoio_writeEsriasciiMap(const string& filename, TInit* pUV, GeoMatrix<double>& gm, long pNumber_novalue){

	Grid2DObject  gridObject;

	if(geotop::common::Variables::UV->U[1]!=geotop::common::Variables::UV->U[2]){
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",geotop::common::Variables::UV->U[2],geotop::common::Variables::UV->U[1]);
		t_error("Fatal error");
	}

	unsigned int ncols = gm.getCols()-1;
	unsigned int nrows = gm.getRows()-1;

	gridObject.grid2D.resize(ncols, nrows, IOUtils::nodata);
    gridObject.llcorner.setXY(geotop::common::Variables::UV->U[4],geotop::common::Variables::UV->U[3], 0);
	gridObject.set(ncols, nrows, geotop::common::Variables::UV->U[1] ,gridObject.llcorner,gridObject.grid2D);

	/* Copies a GeoMatrix to a MeteoIO Grid2DObject */
		for (unsigned int ii = 0; ii < nrows; ii++) {
			for (unsigned int jj = 0; jj < ncols; jj++) {
			//	if (gm(nrows  - ii,jj) == IOUtils::nodata) {
			//		gridObject.grid2D(jj,ii) = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
			//	} else {
					gridObject.grid2D(jj,ii) = gm(nrows -  ii,jj+1);
			//	}
			}
		}

	  io->write2DGrid(gridObject, filename+".asc");
	}

void meteoio_writeEsriasciiVector(const std::string& filenam, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue){

		Grid2DObject  gridObject;

		if(UV->U[1]!=UV->U[2]){
			printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
			t_error("Fatal error");
		}

		gridObject.grid2D.resize(nc, nr, IOUtils::nodata);
		gridObject.llcorner.setXY(UV->U[4],UV->U[3], 0);
	    gridObject.set(nc,nr, UV->U[1] ,gridObject.llcorner,gridObject.grid2D);

		for(long ii=0 ;ii< nr;ii++){
			for(long jj= 0;jj< nc;jj++){
				if (j[nr - ii][jj+1] > 0) {
					if(type==1){
						gridObject.grid2D(jj,ii) = (long)(DTM[j[ii][jj]]);
					}else{
						gridObject.grid2D(jj,ii) = DTM[j[ nr - ii][jj+1]];
					}
				}else {
				//	gridObject.grid2D(jj,ii) =  novalue;
				}
			}
		}
		io->write2DGrid(gridObject, filenam+".asc");
}

/**
 * @brief The following piece of code is to handle snow/rain partition inside the precipitation by correcting
 * the precipitation at each station for the raincorrfact and the snowcorrfact GEOtop part_snow() method 
 * is used to perform the partitions. HNW and TA are read from the MeteoIO MeteoData and then convert TA 
 * to celsius GEOtop part_snow() method is called for partitioning and then rain correction factor and 
 * snow correction factor are added
 *
 * NOTE: THIS FUNCTIONALITY IS IMPLEMENTED IN METEOIO, the appropriate filters need to be configured, e.g.:
 * HNW::filter3    = undercatch_wmo
 * HNW::arg3    = cst 1.3 1
 * where the arguments of the undercatch_wmo filter are "cst {factor for snow} {factor for mixed precipitation}"
 *
 * @param par      Pointer to GEOtop Par class, holding the values for rain and snow correction factors
 * @param meteo    A vector of MeteoData objects (one for each station) at one point in time
 */
void hnw_correction(Par* par, std::vector<mio::MeteoData>& meteo)
{
	double rain=.0, snow=.0, hnw=.0, ta=.0;
	
	for (size_t ii = 0; ii < meteo.size(); ii++) {
		hnw = meteo[ii](MeteoData::HNW);
		
		if (meteo[ii](MeteoData::TA) != IOUtils::nodata) { /* MeteoIO deals with temperature in Kelvin */
			ta = meteo[ii](MeteoData::TA) - 273.15;
		} else {
			ta = meteo[ii](MeteoData::TA);
		}

		/* check for non-zero precipitation, temperature must be available */
		if (hnw > 0 && hnw != IOUtils::nodata && ta != IOUtils::nodata) {
			/* perform pre-processing on HNW using GEOtop part_snow method */
			part_snow(hnw, &rain, &snow, ta, par->T_rain, par->T_snow);

			/* adding rain correction factor and snow correction factor */
			meteo[ii](MeteoData::HNW) = par->raincorrfact * rain + par->snowcorrfact * snow;
		}
	}
}

/**
 * @brief This function performs the 2D interpolations by calling the respective methods within MeteoIO
 * 1) Data is copied from MeteoIO objects to GEOtop structs
 * 2) Interpolation takes place
 * 3) Gridded data copied back to GEOtop DOUBLEMATRIX
 * 4) TA, P and RH values need to be converted as well as nodata values
 *
 * @param par Pointer to the GEOtop Par class, holding geotop.inpts configuration
 * @param currentdate Matlab julian date for which interpolation is desired
 * @param met Pointer to the GEOtop Meteo class 
 * @param wat Pointer to the GEOtop Water class
 */
void meteoio_interpolate(Par* par, double currentdate, Meteo* met, Water* wat) {
	// We need some intermediate storage for storing the interpolated grid by MeteoIO
	Grid2DObject tagrid, rhgrid, pgrid, vwgrid, dwgrid, hnwgrid, cloudwgrid;

	Date d1;
	d1.setMatlabDate(currentdate, geotop::common::Variables::TZ); // GEOtop use matlab offset of julian date 

	try {

		// Intermediate storage for storing data sets for 1 timestep
		std::vector<mio::MeteoData> meteo;

		// Read the meteo data for the given timestep
		io->getMeteoData(d1, meteo);

		//Bypass the internal reading of MeteoData and to performed processing and interpolation 
		//on the data of the given vector meteo
		//hnw_correction(par, meteo);
		//io->add_to_cache(d1, meteo);

		io->getMeteoData(d1, dem, MeteoData::TA, tagrid);
		changeTAgrid(tagrid);

		io->getMeteoData(d1, dem, MeteoData::RH, rhgrid);
		//changeRHgrid(rhgrid);// TODO: check whether GEOtop wants RH from 0 to 100 or from 0 to 1

		io->getMeteoData(d1, dem, MeteoData::P, pgrid);
		changePgrid(pgrid);

		try {
			io->getMeteoData(d1, dem, MeteoData::VW, vwgrid);
		} catch (exception& e) {
		   changeVWgrid(vwgrid, par->Vmin);
		}

		try {
			io->getMeteoData(d1, dem, MeteoData::DW, dwgrid);
		} catch (exception& e) {
		   changeGrid(dwgrid, 90);
		}

		io->getMeteoData(d1, dem, MeteoData::HNW, hnwgrid);

		//io.write2DGrid(vwgrid, "vw_change.asc");

	} catch (exception& e) {
		cerr << "[ERROR] MeteoIO: " << e.what() << endl;
	}

	cout << "[MeteoIO] Start copying Grid to GEOtop format: " << endl;

	// Now copy all that data to the appropriate GEOtop grids
	copyGridToMatrix(tagrid, met->Tgrid);
	copyGridToMatrix(rhgrid, met->RHgrid);
	copyGridToMatrix(pgrid, met->Pgrid);
	copyGridToMatrix(vwgrid, met->Vgrid);
	copyGridToMatrix(dwgrid, met->Vdir);
	copyGridToMatrix(hnwgrid, wat->PrecTot);
}

void replace_grid_values(const DEMObject& dem, const double& value, Grid2DObject& grid) 
{
	grid.set(dem.ncols, dem.nrows, dem.cellsize, dem.llcorner);
	for (unsigned int jj=0; jj<dem.ncols; jj++) {
		for (unsigned int ii=0; ii<dem.nrows; ii++) {
			if (dem(jj,ii) != IOUtils::nodata) {
				grid.grid2D(jj,ii) = value;
			}
		}
	}
}

//void meteoio_interpolate_pointwise(PAR* par, double currentdate,
//		METEO* met, WATER* wat) {
  void meteoio_interpolate_pointwise(Par* par, double currentdate,
			Meteo* met, Water* wat) {

	/*
	 This function performs the Point wise interpolations by calling the respective methods within MeteoIO
	 1) Data is copied from MeteoIO objects to GEOtop structs
	 2) interpolation takes place
	 3) gridded data copied back to GEOtop DOUBLEMATRIX
	 4) TA, P and RH values need to be converted as well as nodata values
	 */

	std::vector<double> resultTa, resultRh, resultP, resultVw, resultDw,
			resultHnw;
	Date d1;
	d1.setMatlabDate(currentdate, geotop::common::Variables::TZ); // GEOtop use matlab offset of julian date

	std::vector<StationData> vecStation;

	try {
		/* Intermediate storage for storing collection of stations data */
		std::vector<std::vector<MeteoData> > vecMeteos;
		/* Intermediate storage for storing data sets for 1 timestep */
		std::vector<MeteoData> meteo;
		/* Read the meteo data for the given timestep */
		size_t numOfStations = io->getMeteoData(d1, meteo);
		/* Allocation for the vectors for storing collection of stations data */
		vecMeteos.insert(vecMeteos.begin(), meteo.size(),std::vector<MeteoData>());

		/*
		 * The following piece of code is to handle snow/rain partition inside the precipitation by correcting
		 * the precipitation at each station for the raincorrfact and the snowcorrfact
		 * GEOtop part_snow() method is used to perform the partitions.
		 * HNW and TA are read from the MeteoIO MeteoData and then convert TA to celsius
		 * GEOtop part_snow() method is called for partitioning and then rain correction factor and snow correction factor are added
		 * Finally push back the modified data to the MeteoIO internal storage
		 */

		double rain, snow, hnw, ta;

		for (size_t i = 0; i < numOfStations  ; i++) {
			hnw = meteo[i](MeteoData::HNW);
			if (meteo[i](MeteoData::TA) != IOUtils::nodata) {
			/* MeteoIO deals with temperature in Kelvin */
				ta = meteo[i](MeteoData::TA) - 273.15;
			} else {
				ta = meteo[i](MeteoData::TA);
			}
			/* check for non-zero precipitation */
			if (hnw > 0 && hnw != IOUtils::nodata && ta != IOUtils::nodata) {
				/* perform pre-processing on HNW using GEOtop part_snow method */
				
				part_snow(hnw, &rain, &snow, ta, par->T_rain, par->T_snow);
				/* adding rain correction factor and snow correction factor */
				meteo[i](MeteoData::HNW) = par->raincorrfact * rain + par->snowcorrfact * snow;
				/* fill the data into the temporary vector of vectors */
				vecMeteos.at(i).push_back(meteo[i]);
				
			}
		}

		/* Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector */
		io->add_to_cache(d1, meteo);

		/* Prepare the points */
		Coords point;
		std::vector<Coords> pointsVec; //create a vector of Coords objects

		double eastX, northY;

		/* Read the X,Y points listed in the GEOtop ListPoint file and prepare a vector to push in the MeteoIO interpolation */

	//	for (int i = 1; i <= par->chkpt->nrh; i++) {
		for (size_t i = 1; i < par->chkpt.getRows(); i++) {
		//	eastX = par->chkpt->co[i][ptX];
			eastX = par->chkpt[i][ptX];
		//	northY = par->chkpt->co[i][ptY];
			northY = par->chkpt[i][ptY];
			point.setXY(eastX, northY, IOUtils::nodata);
			pointsVec.push_back(point);
		}

		std::cout << "\n[MeteoIO] Time to interpolate Point wise: "
					<< d1.toString(Date::ISO) << std::endl;

		/*Push the coordinate list in the MeteoIO interpolate method tom perform the interpolation on them */

		io->interpolate(d1, dem, MeteoData::TA, pointsVec, resultTa);
		/* change TA values */
		for (size_t i = 0; i < resultTa.size(); i++) {
			if (resultTa[i] != IOUtils::nodata) {
				resultTa[i] -= 273.15; //back to celsius
			}
		}

		io->interpolate(d1, dem, MeteoData::RH, pointsVec, resultRh);

		io->interpolate(d1, dem, MeteoData::P, pointsVec, resultP);
		/* change P values */
		for (size_t i = 0; i < resultP.size(); i++) {
			if (resultP[i] != IOUtils::nodata) {
				resultP[i] /= 100.0;
			}
		}

		try {
			io->interpolate(d1, dem, MeteoData::VW, pointsVec, resultVw);
		} catch (std::exception& e) {
			/* change VW values, if the measured value is less than the user given vmin value then set the vmin */
			for (size_t i = 0; i < resultVw.size(); i++) {
				if (resultVw[i] != IOUtils::nodata && resultVw[i] < par->Vmin) {
					resultVw[i] = par->Vmin; //set GEOtop vw min value
				}
			}
		}
		try{
			io->interpolate(d1, dem, MeteoData::DW, pointsVec, resultDw);
		}catch (std::exception& e){
			/* change DW values to a user given  value */
			for (size_t i = 0; i < resultDw.size(); i++) {
				if (resultDw[i] != IOUtils::nodata) {
					resultDw[i] = 90;
				}
			}
		}

		io->interpolate(d1, dem, MeteoData::HNW, pointsVec, resultHnw);


	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	cout << "[MeteoIO] Start copying point data to GEOtop format: " << endl;
	/* Now copy all that MeteoIO interpolated data to the appropriate GEOtop grids */

	copyGridToMatrixPointWise(resultTa, met->Tgrid);
	copyGridToMatrixPointWise(resultRh, met->RHgrid);
	copyGridToMatrixPointWise(resultP, met->Pgrid);
	copyGridToMatrixPointWise(resultDw, met->Vdir);
	copyGridToMatrixPointWise(resultVw, met->Vgrid);
	copyGridToMatrixPointWise(resultHnw, wat->PrecTot);
}

void copyInterpMeteoData(double *out, std::vector<MeteoData>& meteoin){

	for(size_t i=0;i<meteoin.size();i++)
	{
//	out[0]=151020111700.000000;   // timestamp
//	out[1]=734791.666667;  		  // timestamp

	out[2]=checkNOvalue(meteoin[i](MeteoData::HNW));
	out[3]=checkNOvalue(IOUtils::nodata);
	out[4]=checkNOvalue(meteoin[i](MeteoData::VW));
	out[5]=checkNOvalue(meteoin[i](MeteoData::DW));
	if(meteoin[i](MeteoData::VW)!=IOUtils::nodata && meteoin[i](MeteoData::DW) !=IOUtils::nodata){
		out[6]=-meteoin[i](MeteoData::VW) * sin(meteoin[i](MeteoData::DW) * GTConst::Pi / 180.);
		out[7]=-meteoin[i](MeteoData::VW)* cos(meteoin[i](MeteoData::DW) * GTConst::Pi / 180.);
	} else{
		out[6]=geotop::input::gDoubleNoValue;
		out[7]=geotop::input::gDoubleNoValue;
	}

        out[8]=meteoin[i](MeteoData::RH)==IOUtils::nodata? geotop::input::gDoubleNoValue : meteoin[i](MeteoData::RH)*100;
	out[9]= meteoin[i](MeteoData::TA)==IOUtils::nodata? geotop::input::gDoubleNoValue : meteoin[i](MeteoData::TA)-273.15;

	if(meteoin[i](MeteoData::TA)!=IOUtils::nodata && meteoin[i](MeteoData::RH)!=IOUtils::nodata && meteoin[i](MeteoData::P)!=IOUtils::nodata)
	{
		out[10]=tDew(out[9], out[8], meteoin[i](MeteoData::P) /100.0);// see Tdew(double T, double RH, double P)
	}else{
		out[10]=geotop::input::gDoubleNoValue;
	}

	out[11]=checkNOvalue(meteoin[i](MeteoData::ISWR));
	out[12]=geotop::input::gDoubleNoValue;// SWb (beam component)
	out[13]=geotop::input::gDoubleNoValue;// SWd (diffuse component)
	out[14]=geotop::input::gDoubleNoValue;// cloud transmissivity
	out[15]=geotop::input::gDoubleNoValue;// cloud factor
	out[16]=geotop::input::gDoubleNoValue;// LWin
	out[17]=geotop::input::gDoubleNoValue;// SWnet

//	cout<<"Station "<<i+1<<endl;
//	cout<<"meteo[2] " <<out[2]<<endl;
//	cout<<"meteo[4] " <<out[4]<<endl;
//	cout<<"meteo[5] " <<out[5]<<endl;
//	cout<<"meteo[6] " <<out[6]<<endl;
//	cout<<"meteo[7] " <<out[7]<<endl;
//	cout<<"meteo[8] " <<out[8]<<endl;
//	cout<<"meteo[9] " <<out[9]<<endl;
//	cout<<"meteo[10] "<<out[10]<<endl;
//	cout<<"meteo[11] "<<out[11]<<endl;
//	cout<<"P " <<checkNOvalue(meteoin[i](MeteoData::P))<<endl;
	}

}

double checkNOvalue(double var){
	return var==IOUtils::nodata? geotop::input::gDoubleNoValue : var;
}

double tDew(double T, double RH, double P){
	/*
	satVapPressure: vapour pressure [mbar]
	P: atmospheric pressure [mbar]
	T: air temperature in [C]
	*/
	double satVapPressure, tFromSatVapPressure, A, b, c;

	A=6.1121*(1.0007+3.46E-6*P);
	b=17.502;
	c=240.97;

	//calculate saturated vapour pressure
	satVapPressure=A*exp(b*T/(c+T));

	//find min and multiply it with P
	satVapPressure *=(RH <1.0 ? RH : 1.0);

	//calculate temperature from saturated water vapour
	tFromSatVapPressure = c*log(satVapPressure/A)/(b-log(satVapPressure/A));

	return tFromSatVapPressure;
}


/*void copyGridToMatrix(Grid2DObject& gridObject, DOUBLEMATRIX* myGrid) {

	 This function copy MeteoIO Grid2DObject to the GEOtop structs

	for (size_t ii = 0; ii < gridObject.nrows; ii++) {
		for (size_t jj = 0; jj < gridObject.ncols; jj++) {
			if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii) == IOUtils::nodata) {
				myGrid->co[ii + 1][jj + 1] = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
			} else {
				myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj, gridObject.nrows - 1 - ii);
			}
			//std::cout<<"myGrid->co["<<ii<<"]["<<jj<<"]"<<myGrid->co[ii + 1][jj + 1] << std::endl;
		}
	}
}*/

//overloaded method- noori
void copyGridToMatrix(Grid2DObject& gridObject, GeoMatrix<double>& myGrid) {

	/* This function copy MeteoIO Grid2DObject to the GEOtop structs */

	for (size_t ii = 0; ii < gridObject.nrows; ii++) {
		for (size_t jj = 0; jj < gridObject.ncols; jj++) {
			if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii) == IOUtils::nodata) {
				myGrid[ii + 1][jj + 1] = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
			} else {
				myGrid[ii + 1][jj + 1] = gridObject.grid2D(jj, gridObject.nrows - 1 - ii);
			}
			//std::cout<<"myGrid->co["<<ii<<"]["<<jj<<"]"<<myGrid->co[ii + 1][jj + 1] << std::endl;
		}
	}
}


void copyGridToMatrixPointWise(std::vector<double>& pointValues,
		GeoMatrix<double>& myGrid) {

	/*
	 * This function copy MeteoIO Point vector to the GEOtop structs
	 * Notice that, co[1][1] is the first index to accessed, because in GEOtop there is an index offset
	 * therefor we have to start with [1][1]
	 */

	for (unsigned int i = 0; i < pointValues.size(); i++) {
		if (pointValues[i] != IOUtils::nodata) {
			myGrid[1][i + 1] = pointValues[i]; //co[1][1] is the first index accessed
		} else {
			myGrid[1][i + 1] = geotop::input::gDoubleNoValue; //co[1][1] is the first index accessed
		}
	}
}


void changeRHgrid(Grid2DObject& g2d) {

	/* This function changes the RH according to to GEOtop */

	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) *= 100.0;
		}
	}
}

void changeTAgrid(Grid2DObject& g2d) {

	/* This function changes the TA, kelvin to celsius according to to GEOtop */

	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) -= 273.15; //back to celsius
		}
	}
}

void changePgrid(Grid2DObject& g2d) {

	/* This function changes the P according to to GEOtop */

	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) /= 100.0;
		}
	}
}

void changeVWgrid(Grid2DObject& g2d, double vwMin) {

	/* This function changes the VW according to the GEOtop given min value */

	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata
					&& g2d.grid2D(ii, jj) < vwMin)
				g2d.grid2D(ii, jj) = vwMin;
		}
	}
}

void changeGrid(Grid2DObject& g2d, const double val) {
	/* This function changes a grid according to given value */
	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) = val;
		}
	}
}

//double ***meteoio_readMeteoData(long*** column,
//		METEO_STATIONS *stations, long nrOfVariables, PAR *par,
//		TIMES *times) {

double ***meteoio_readMeteoData(long*** column,
		MeteoStations *stations, long nrOfVariables, Par *par,
		Times *times) {


	long ncols = nrOfVariables; //the total number of meteo variables used in GEOtop (should stay fixed)

	//Date d1 holds the beginning of this simulation, d2 the end date of the simulation
	//d1=times->time;
	//d2=times->time+par->Dt;

//	Date d1 = par->init_date->co[geotop::common::Variables::i_sim];
	Date d1 = par->init_date[geotop::common::Variables::i_sim];
//	Date d2 = par->end_date->co[geotop::common::Variables::i_sim];
	Date d2 = par->end_date[geotop::common::Variables::i_sim];

	//	Date d1((int)par->year0, 1, 1, 0, 0);
	//	d1 += par->JD0;
	//	d1 += times->time/86400; //times->time is in seconds, conversion to julian by devision
	//
	//	Date d2((int)par->year0, 1, 1, 0, 0);
	//	d2 += par->JD0;
	//	d2 += times->TH/24;      //the end of the simulation

	//Construction a BufferIOHandler and reading the meteo data through meteoio as configured in io.ini
	Config cfg("io.ini");
	//IOManager io(cfg);

	std::vector<std::vector<MeteoData> > vecMeteo; //the dimension of this vector will be nrOfStations
	std::vector<std::vector<StationData> > vecStation; //the dimension of this vector will be nrOfStations

	cout << "[I] MeteoIO: Fetching all data between " << d1.toString(Date::ISO)
			<< " and " << d2.toString(Date::ISO) << endl;

	/*
	 * In the following section the BufferedIOHandler is used to request data for descrete times
	 * MeteoIO will deal with accumulation, linear interpolation, filtering and cleaning of the data
	 * During the loop a discrete amount of time will be added to the loop variable (one hour)
	 */
	for (Date currentDate = d1; currentDate <= d2;
			currentDate += Date(1.0 / 24.0)) {
		cout << "[I] MeteoIO: Getting data for "
				<< currentDate.toString(Date::ISO) << endl;
		std::vector<MeteoData> vecM; //dimension: nrOfStations
		std::vector<StationData> vecS; //dimension: nrOfStations

		io->getMeteoData(currentDate, vecM); //reading meteo data and meta datafor currentDate

		//the data needs to be appended to the collection of all read meteo and meta data:
		for (unsigned int jj = 0; jj < vecM.size(); jj++) {
			if (currentDate == d1) {
				vecMeteo.push_back(std::vector<MeteoData>());
				vecStation.push_back(std::vector<StationData>());
			}
			vecMeteo.at(jj).push_back(vecM.at(jj)); //append meteo data to vector
			vecStation.at(jj).push_back(vecS.at(jj)); //append meta data to vector
		}
	}
	cout << "[I] MeteoIO: " << "Finished getting meteo and station data"
			<< endl;
	//print out all data if configured by user
	try {
		string tmp;
		cfg.getValue("METEO", "OUTPUT", tmp);
		io->writeMeteoData(vecMeteo);
	} catch (std::exception& e) {
	} //Do nothing if not configured or error happens

	//Deal with meta data, that is allocate met->st and fill with data
	std::vector<StationData> vecMyStation;
	io->getStationData(d1, vecMyStation);

	initializeMetaData(vecMyStation, d1,  par, stations);
	cout << "[I] MeteoIO: Amount of stations: " << vecMeteo.size() << endl;

	//resize the column matrix - it will hold information about what column in the data matrix holds what parameter
	*column = (long**) realloc(*column, vecStation.size() * sizeof(long*));

	double ***data = (double ***) malloc(vecMeteo.size() * sizeof(double**));
	short novalueend = 1;

	for (unsigned int jj = 0; jj < vecStation.size(); jj++) { //for each station
		//{Iprec, WindS, WindDir, RelHum, AirT, AirP, SWglobal, SWdirect, SWdiffuse, TauCloud, Cloud, LWin, SWnet, Tsup}
		//iPt ,iWs ,iWdir , iWsx , iWsy, iRh ,iT ,iTdew ,iPs,iSW ,iSWb ,iSWd,itauC ,iC,iLWi,iSWn
		(*column)[jj] = (long*) malloc((ncols + 1) * sizeof(long));
		//(*column)[jj][ncols] = end_vector_long;
		for (int ff = 0; ff < nrOfVariables; ff++) {
			//(*column)[jj][ff] = ff;
			(*column)[jj][ff] = -1;
		}

		unsigned int ii = 0;
		while ((*column)[jj][ii] != 999999) { //the geotop way of ending a vector
			//cout << "Station " << jj << "  " << ii << ": " << (*column)[jj][ii] << endl;
			ii++;
		}

		data[jj] = (double **) malloc(vecMeteo[jj].size() * sizeof(double*));

		unsigned int ll = 0;
		for (ll = 0; ll < vecMeteo[jj].size(); ll++) {
			data[jj][ll] = (double *) malloc((ncols + 1) * sizeof(double));

			//Put data in the correct cell

			data[jj][ll][0] = vecMeteo[jj][ll].HNW;
			data[jj][ll][1] = vecMeteo[jj][ll].VW;
			data[jj][ll][2] = vecMeteo[jj][ll].DW;
			data[jj][ll][3] = vecMeteo[jj][ll].RH * 100.00; //MeteoIO deals with RH in values [0;1]
			data[jj][ll][4] = vecMeteo[jj][ll].TA - 273.15; //MeteoIO deals with temperature in Kelvin
			data[jj][ll][5] = vecMeteo[jj][ll].P;
			data[jj][ll][6] = vecMeteo[jj][ll].ISWR;
			data[jj][ll][7] = geotop::input::gDoubleNoValue;
			data[jj][ll][8] = geotop::input::gDoubleNoValue;
			data[jj][ll][9] = geotop::input::gDoubleNoValue;
			data[jj][ll][10] = geotop::input::gDoubleNoValue;
			data[jj][ll][11] = geotop::input::gDoubleNoValue;
			data[jj][ll][12] = geotop::input::gDoubleNoValue;
			data[jj][ll][13] = geotop::input::gDoubleNoValue;
			//data[jj][ll][ncols] = end_vector;

			for (int gg = 0; gg < nrOfVariables; gg++) {
				if (data[jj][ll][gg] == IOUtils::nodata) {
					data[jj][ll][gg] = geotop::input::gDoubleNoValue;
				} else if (data[jj][ll][gg] != geotop::input::gDoubleNoValue) {
					(*column)[jj][gg] = gg; //HACK!
				}
			}

			//Comment this in if you dont want to see the data that was retrieved
			cout << "MeteoData [" << ll + 1 << "/" << vecMeteo[jj].size()
					<< "]:" << endl;
			//	std::cout << vecMeteo[jj][ll].toString() << std::endl;

		}

		for (int ff = 1; ff <= ncols; ff++) {
			if (ll > 0)
				if (data[jj][ll - 1][ff] != geotop::input::gDoubleNoValue)
					novalueend = 0;
		}

		//		if(novalueend==0){
		//			data[jj]=(double **)realloc(data[jj],(vecMeteo[jj].size()+1)*sizeof(double*));
		//			data[jj][vecMeteo[jj].size()] = (double *)malloc(sizeof(double));
		//			//data[jj][vecMeteo[jj].size()][0]=end_vector;
		//		} else {
		//			if (vecMeteo[jj].size()>0)
		//				data[jj][vecMeteo[jj].size()-1][0]=end_vector;
		//		}
	}

	cout << "[I] MeteoIO NOVAL used          : " << geotop::input::gDoubleNoValue << endl;
	cout << "[I] MeteoIO #of meteo parameters: " << nrOfVariables << endl;

	//Testing access to the whole tensor
	for (size_t ii = 0; ii < vecStation.size(); ii++) {
		double d = 0.0;
		for (size_t ll = 0; ll < vecMeteo[ii].size(); ll++) {
			for (int kk = 0; kk < ncols; kk++) {
				//printf("%f ", data[ii][ll][kk]);
				d = data[ii][ll][kk];
			}
			//printf("\n");
		}
	}

	return data; //return pointer to heap allocated memory
}

//void initializeMetaData(const std::vector<StationData>& vecStation,
//		const Date& startDate,  PAR *par,
//		METEO_STATIONS *stations) {

void initializeMetaData(const std::vector<StationData>& vecStation,
		const Date& startDate,  Par *par,
		MeteoStations *stations) {
	unsigned int nrOfStations = vecStation.size();

	//Initialize the station data: set beginning of data
	int year, month, day, hour, minute;
	startDate.getDate(year, month, day, hour, minute);
	//	Date JD_start = startDate - Date((int)par->year0, 1, 1, 0, 0); // this is the actual elapsed simulation time

	//init station struct met->st
//	stations->E = new_doublevector(nrOfStations); // East coordinate [m] of the meteo station
	stations->E.resize(nrOfStations+1); // East coordinate [m] of the meteo station
//	stations->N = new_doublevector(nrOfStations); // North coordinate [m] of the meteo station
	stations->N.resize(nrOfStations+1); // North coordinate [m] of the meteo station
//	stations->lat = new_doublevector(nrOfStations); // Latitude [rad] of the meteo station
	stations->lat.resize(nrOfStations+1); // Latitude [rad] of the meteo station
//	stations->lon = new_doublevector(nrOfStations); // Longitude [rad] of the meteo station
	stations->lon.resize(nrOfStations+1); // Longitude [rad] of the meteo station
//	stations->Z = new_doublevector(nrOfStations); // Elevation [m] of the meteo station
	stations->Z.resize(nrOfStations+1); // Elevation [m] of the meteo station
//	stations->sky = new_doublevector(nrOfStations); // Sky-view-factor [-] of the meteo station
	stations->sky.resize(nrOfStations+1); // Sky-view-factor [-] of the meteo station
//	stations->ST = new_doublevector(nrOfStations); // Standard time minus UTM [hours] of the meteo station
	stations->ST.resize(nrOfStations+1); // Standard time minus UTM [hours] of the meteo station
//	stations->Vheight = new_doublevector(nrOfStations); // Wind velocity measurement height [m] (a.g.l.)
	stations->Vheight.resize(nrOfStations+1); // Wind velocity measurement height [m] (a.g.l.)
//	stations->Theight = new_doublevector(nrOfStations); // Air temperature measurement height [m] (a.g.l.)
	stations->Theight.resize(nrOfStations+1); // Air temperature measurement height [m] (a.g.l.)

	//	stations->JD0=new_doublevector(nrOfStations);// Decimal Julian Day of the first data
	//	stations->Y0=new_longvector(nrOfStations);	// Year of the first data
	//	stations->Dt=new_doublevector(nrOfStations);	// Dt of sampling of the data [sec]
	//	stations->offset=new_longvector(nrOfStations);// offset column

	for (unsigned int ii = 1; ii <= nrOfStations; ii++) { //HACK
		std::cout << "[I] MeteoIO station " << ii << ":\n" << vecStation[ii - 1].toString()
				<< std::endl;

	//	stations->E->co[ii] = vecStation[ii - 1].position.getEasting();
		stations->E[ii] = vecStation[ii - 1].position.getEasting();
	//	stations->N->co[ii] = vecStation[ii - 1].position.getNorthing();
		stations->N[ii] = vecStation[ii - 1].position.getNorthing();
	//	stations->lat->co[ii] = vecStation[ii - 1].position.getLat() * PI/ 180.0; // from deg to [rad]
		stations->lat[ii] = vecStation[ii - 1].position.getLat() * PI/ 180.0; // from deg to [rad]
	//	stations->lon->co[ii] = vecStation[ii - 1].position.getLon() * PI/ 180.0; // from deg to [rad]
		stations->lon[ii] = vecStation[ii - 1].position.getLon() * PI/ 180.0; // from deg to [rad]
	//	stations->Z->co[ii] = vecStation[ii - 1].position.getAltitude();
		stations->Z[ii] = vecStation[ii - 1].position.getAltitude();
	//	stations->sky->co[ii] = geotop::input::gDoubleNoValue;
		stations->sky[ii] = geotop::input::gDoubleNoValue;
	//	stations->ST->co[ii] = 0;
		stations->ST[ii] = 0;
	//	stations->Vheight->co[ii] = 8.0;
		stations->Vheight[ii] = 8.0;
	//	stations->Theight->co[ii] = 8.0;
		stations->Theight[ii] = 8.0;

	}
}
/**
 * @brief Function that checks in a vector of MeteoData if a valid ISWR value is present
 *
 * @param vec_meteo vector of MeteoData
 * @param first_check true if this is the first check, a special warning message will be printed
 * @param A The AllData reference
 *
 * @return Return true if at least one station measured ISWR 
 */
bool iswr_present(const std::vector<mio::MeteoData>& vec_meteo, const bool& first_check, AllData *A) 
{
	A->M->nstcloud   = 0; // first meteo station ID (1...n) to use for the cloudiness (ISWR)
	A->M->numstcloud = 0; // counter of meteo stations containing cloud info

	for (size_t ii=0; ii<vec_meteo.size(); ii++) {
		if (vec_meteo[ii](MeteoData::ISWR) != IOUtils::nodata) { //ISWR is measured
			if(A->M->nstcloud==0) A->M->nstcloud = ii; // records the first station that has cloud info

			A->M->numstcloud++;// increase the counter of meteo stations measuring cloud

			//	A->M->st->flag_SW_meteoST->co[ii] = 1;
			A->M->st->flag_SW_meteoST[ii] = 1; //HACK: this index should be checked
		}
	}
	
	bool iswr_found = (A->M->numstcloud > 0);

	if (first_check) {
		if (!iswr_found) {
			printf("WARNING: NO cloudiness measurements available\n");
			//fprintf(flog,"WARNING: NO cloudiness measurements available\n"); //HACK

			A->M->nstcloud = 1;
		} else {
			printf("Cloudiness measurements from %ld stations: will be used data from station %ld\n",A->M->numstcloud,A->M->nstcloud);
			//fprintf(flog,"Cloudiness measurements from station %ld\n",met->nstcloud); //HACK
		}
	}

	return iswr_found;
}

//void meteoio_interpolate_cloudiness(T_INIT* pUV, PAR* par,
//		double currentdate, DOUBLEMATRIX* tau_cloud_grid,
//		DOUBLEVECTOR* tau_cloud_vec) {

void meteoio_interpolate_cloudiness(TInit* pUV, Par* par,
		double currentdate, GeoMatrix<double>& tau_cloud_grid,
		GeoVector<double>& tau_cloud_vec) {
	/*
	 This function performs the 2D interpolations of current cloudiness by calling the respective methods within MeteoIO
	 1) Data is copied from MeteoIO objects to GEOtop structs
	 2) interpolation takes place
	 4) gridded data copied back to GEOtop DOUBLEMATRIX inside the tau_cloud_grid
	 */

	Grid2DObject cloudwgrid;

	Date d1;
	d1.setMatlabDate(currentdate, geotop::common::Variables::TZ); // GEOtop use matlab offset of julian date

	std::cout << "\n[MeteoIO] Time to interpolate : cloudiness "
			<< d1.toString(Date::ISO) << std::endl;

	try {
		std::vector<std::vector<MeteoData> > vecMeteos;
		std::vector<MeteoData> meteo; // Intermediate storage for storing data sets for 1 timestep
		int numOfStations = io->getMeteoData(d1, meteo); // Read the meteo data for the given timestep
		vecMeteos.insert(vecMeteos.begin(), meteo.size(), std::vector<MeteoData>()); // Allocation for the vectors

		for (int i = 0; i < numOfStations; i++) {
		//	meteo[i](MeteoData::RSWR) = tau_cloud_vec->co[i];
			meteo[i](MeteoData::RSWR) = tau_cloud_vec[i];
			vecMeteos.at(i).push_back(meteo[i]); // fill the data into the vector of vectors
		}

		// Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector
		io->push_meteo_data(IOManager::filtered, d1, d1, vecMeteos);

		// if point_sim == 1 then point-wise simulation otherwise distributed simulation

		if (par->point_sim == 1) {

			//prepare the points
			Coords point;
			std::vector<Coords> pointsVec; //create a vector of Coords objects
			std::vector<double> resultCloud;
			double eastX, northY;

		//	for (int i = 1; i <= par->chkpt->nrh; i++) {
			for (size_t i = 1; i < par->chkpt.getRows(); i++) {
			//	eastX = par->chkpt->co[i][ptX];
				eastX = par->chkpt[i][ptX];
			//	northY = par->chkpt->co[i][ptY];
				northY = par->chkpt[i][ptY];
				point.setXY(eastX, northY, IOUtils::nodata);
				pointsVec.push_back(point);
			}
            /* Interpolate point wise */
			try {
				io->interpolate(d1, dem, MeteoData::RSWR, pointsVec, resultCloud);
			} catch (std::exception& e) {
				for (size_t i = 0; i < resultCloud.size(); i++) {
					if (resultCloud[i] != IOUtils::nodata) {
						resultCloud[i] = 0.5;
					}
				}
			}
			/* Now copy all that data to the appropriate Array */
			copyGridToMatrixPointWise(resultCloud, tau_cloud_grid);


		} else {
			/*Interpolate 2D grid*/
			try {
				io->getMeteoData(d1, dem, MeteoData::RSWR, cloudwgrid);
			} catch (std::exception& e) {
			   changeGrid(cloudwgrid, 0.5);
			}
			/* Now copy all that data to the appropriate grids */
			copyGridToMatrix(cloudwgrid, tau_cloud_grid);
		}

	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}
}


		// ---------   read the GEOtop lapseRate and set it to MeteoIO config -----

		//		double lrta = -LapseRateTair * 1E-3; // GEOtop default lapse rate of air temp.
		//		double lspr = -LapseRatePrec * 1E-3; // GEOtop default lapse rate of precipitation
		//
		//		if (met->LRv[ilsTa] != LapseRateTair) {
		//			// The user has given a lapse rate: use this
		//			lrta = -met->LRv[ilsTa] * 1E-3;
		//			stringstream lapserateTA;	// Create a stringstream
		//			lapserateTA << lrta;		// Add number to the stream
		//			string s = " soft";
		//			cfg.addKey("TA::idw_lapse", "Interpolations2D", lapserateTA.str()
		//					+ s);
		//		}
		//
		//		if (met->LRv[ilsPrec] != LapseRatePrec) {
		//			// The user has given a lapse rate: use this
		//			lspr = -met->LRv[ilsPrec] * 1E-3;
		//			stringstream lapserateHNW;	// Create a stringstream
		//			lapserateHNW << lspr;		// Add number to the stream
		//			string s = " soft";
		//			cfg.addKey("HNW::idw_lapse", "Interpolations2D", lapserateHNW.str()
		//					+ s);
		//		}
		//
		//		IOManager io(cfg);
		// -------------End of reconfig io.ini----------------------

