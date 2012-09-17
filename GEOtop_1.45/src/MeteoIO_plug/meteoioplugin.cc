#ifdef USE_METEOIO
#include "meteoioplugin.h"
#include "../geotop/times.h"
#include <sstream>
#include <stdlib.h>

using namespace std;
using namespace mio;
extern long i_sim;
const double TZ = 1.;
extern long number_novalue;

DEMObject dem;
Config cfg("io_it.ini");
IOManager io(cfg);

extern "C" void meteoio_init() {

	/*
	 *  This function performs the initialization by calling the respective methods within MeteoIO
	 *  1) Reading the DEM (necessary for several spatial interpolations algorithms)
	 */

	try {
		//read the DEM;
		io.readDEM(dem);
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}
}

extern "C" DOUBLEMATRIX *meteoio_readDEM(T_INIT** UVREF) {
	DOUBLEMATRIX *myDEM = NULL;

	T_INIT* UV = (T_INIT *) malloc(sizeof(T_INIT));

	try {

		//	DEMObject dem;
		IOHandler iohandler(cfg);
		iohandler.readDEM(dem);

		myDEM = new_doublematrix(dem.nrows, dem.ncols);
		UV->V = new_doublevector(2);
		UV->U = new_doublevector(4);

		UV->V->co[1] = -1.0;
		UV->V->co[2] = number_novalue; //GEOtop nodata -9999.0

		UV->U->co[1] = dem.cellsize;
		UV->U->co[2] = dem.cellsize;
		UV->U->co[3] = dem.llcorner.getNorthing();
		UV->U->co[4] = dem.llcorner.getEasting();

		for (unsigned int ii = 0; ii < dem.nrows; ii++) {
			for (unsigned int jj = 0; jj < dem.ncols; jj++) {
				if (dem.grid2D(jj, dem.nrows - 1 - ii) == IOUtils::nodata) {
					myDEM->co[ii + 1][jj + 1] = UV->V->co[2];
				} else {
					myDEM->co[ii + 1][jj + 1] = dem.grid2D(jj,
							dem.nrows - 1 - ii);
				}
			}
		}
		std::cout << "Read DEM:" << dem.nrows << "(rows) " << dem.ncols
				<< "(cols)" << std::endl;
		free((*UVREF));
		(*UVREF) = UV;
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	return myDEM;
}

extern "C" DOUBLEMATRIX *meteoio_read2DGrid(T_INIT* UV, char* _filename) {
	DOUBLEMATRIX *myGrid = NULL;

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

		if (UV->U->co[1] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (UV->U->co[2] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (UV->U->co[3] != gridObject.llcorner.getNorthing())
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (UV->U->co[4] != gridObject.llcorner.getEasting())
			throw IOException("Inconsistencies between 2D Grids read", AT);

		for (unsigned int ii = 0; ii < gridObject.nrows; ii++) {
			for (unsigned int jj = 0; jj < gridObject.ncols; jj++) {
				if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii)
						== IOUtils::nodata) {
					myGrid->co[ii + 1][jj + 1] = UV->V->co[2];
				} else {
					myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj,
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

	return myGrid;
}

extern "C" void meteoio_interpolate(PAR* par, double currentdate, METEO* met,
		WATER* wat) {
	/* 
	 This function performs the 2D interpolations by calling the respective methods within MeteoIO
	 1) Data is copied from MeteoIO objects to GEOtop structs
	 2) Interpolation takes place
	 3) Gridded data copied back to GEOtop DOUBLEMATRIX
	 4) TA, P and RH values need to be converted as well as nodata values
	 */


	/* We need some intermediate storage for storing the interpolated grid by MeteoIO */
	Grid2DObject tagrid, rhgrid, pgrid, vwgrid, dwgrid, hnwgrid, cloudwgrid;

	Date d1;
	/* GEOtop use matlab offset of julian date */
	d1.setMatlabDate(currentdate, TZ);

	try {

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

		/* Intermediate storage for storing collection of stations data */
		std::vector<std::vector<MeteoData> > vecMeteos;
		/* Intermediate storage for storing data sets for 1 timestep */
		std::vector<MeteoData> meteo;
		/* Read the meteo data for the given timestep */
		size_t numOfStations = io.getMeteoData(d1, meteo);
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

		// check if all station have NA for HNW
		bool HNW_allNA=true;
		for (size_t i = 0; i < numOfStations; i++) {
			if (meteo[i](MeteoData::HNW) != IOUtils::nodata) {
				HNW_allNA=false;
				break;
			}
		}
		if (HNW_allNA==true) {
			for (size_t i = 0; i < numOfStations; i++) {
				if (meteo[i](MeteoData::HNW) == IOUtils::nodata) {
					meteo[i](MeteoData::HNW) = 0.;
					break;
				}
			}
		}

		for (size_t i = 0; i < numOfStations; i++) {
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
				meteo[i](MeteoData::HNW) = par->raincorrfact * rain
						+ par->snowcorrfact * snow;
				/* fill the data into the temporary vector of vectors */
				vecMeteos.at(i).push_back(meteo[i]);
			}
		}

		/* Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector */
		io.add_to_cache(d1, meteo);

		/* Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector */
		//io.push_meteo_data(IOManager::filtered, d1, d1, vecMeteos);

		// PRECIPITATION
		// TODO: correct the precipitation at each station for the raincorrfact and the snowcorrfact
		//		   double prec, rain, snow;
		//		   for(n=1;n<=met->st->Z->nh;n++){
		//					if((long)met->var[n-1][Pcode]!=number_novalue && (long)met->var[n-1][Pcode]!=number_absent){// check if exists prec. value
		//					prec = met->var[n-1][Pcode];// precipitation of the meteo station n
		//					part_snow(prec, &rain, &snow, met->var[n-1][Tcode], par->T_rain, par->T_snow);
		//					met->var[n-1][Pcode] = par->raincorrfact * rain + par->snowcorrfact * snow;
		//					//hnwgrid.grid2D(r, c) = par->raincorrfact * rain + par->snowcorrfact * snow;
		//					}
		//				}
		//		    part_snow(prec, &rain, &snow, tagrid.grid2D(r, c), Train, Tsnow);
		//
		io.interpolate(d1, dem, MeteoData::TA, tagrid);
		changeTAgrid(tagrid);//convert values accordingly, necessary for TA , RH and P
		io.interpolate(d1, dem, MeteoData::RH, rhgrid);
		io.interpolate(d1, dem, MeteoData::P, pgrid);
		changePgrid(pgrid);//convert values accordingly, necessary for TA , RH and P
		io.interpolate(d1, dem, MeteoData::VW, vwgrid);
//		io.write2DGrid(vwgrid, "vwgrid"+d1.toString(Date::NUM)+".asc");
		changeVWgrid(vwgrid, par->Vmin);//convert values accordingly, necessary for TA , RH and P
		io.interpolate(d1, dem, MeteoData::DW, dwgrid);
		io.interpolate(d1, dem, MeteoData::HNW, hnwgrid);


		//convert values accordingly, necessary for TA , RH and P
		//changeRHgrid(rhgrid);// TODO: check whether GEOtop wants RH from 0 to 100 or from 0 to 1
		//changeTAgrid(tagrid);
		//Date camparedate1();
		//if ((d1 >= comparedate1) && (d1 <= comparedate2)){
			//io.write2DGrid(hnwgrid, "hnw"+d1.toString(Date::NUM)+".asc");
			//}
		//io.write2DGrid(vwgrid, "vw_change.asc");

	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	std::cout << "[MeteoIO] Start copying Grid to GEOtop format: " << std::endl;

	// Now copy all that data to the appropriate grids
	copyGridToMatrix(tagrid, met->Tgrid);
	copyGridToMatrix(rhgrid, met->RHgrid);
	copyGridToMatrix(pgrid, met->Pgrid);
	copyGridToMatrix(vwgrid, met->Vgrid);
	copyGridToMatrix(dwgrid, met->Vdir);
	copyGridToMatrix(hnwgrid, wat->PrecTot);
}

extern "C" void meteoio_interpolate_pointwise(PAR* par, double currentdate,
		METEO* met, WATER* wat) {
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
	d1.setMatlabDate(currentdate, TZ); // GEOtop use matlab offset of julian date

	std::vector<StationData> vecStation;

	try {
		/* Intermediate storage for storing collection of stations data */
		std::vector<std::vector<MeteoData> > vecMeteos;
		/* Intermediate storage for storing data sets for 1 timestep */
		std::vector<MeteoData> meteo;
		/* Read the meteo data for the given timestep */
		size_t numOfStations = io.getMeteoData(d1, meteo);
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

		// matteo: added
//		for (size_t i = 0; i < numOfStations  ; i++) {
//			double Iprec =meteo[i](MeteoData::HNW);
//			double airT =meteo[i](MeteoData::TA)- 273.15;;
//			if (Iprec > 0 && airT < 0 && Iprec != IOUtils::nodata && airT != IOUtils::nodata) {
//				std::cout << "Tset1: Before correction St."<<i+1<<" HNW: " << Iprec <<" TA:"<<airT<<endl;
//			}
//		}
//		stop_execution();
		// matteo added


		// check if all station have NA for HNW
		bool HNW_allNA=true;
		for (size_t i = 0; i < numOfStations; i++) {
			if (meteo[i](MeteoData::HNW) != IOUtils::nodata) {
				HNW_allNA=false;
				break;
			}
		}
		if (HNW_allNA==true) {
			for (size_t i = 0; i < numOfStations; i++) {
				if (meteo[i](MeteoData::HNW) == IOUtils::nodata) {
					meteo[i](MeteoData::HNW) = 0.;
					break;
				}
			}
		}


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
						meteo[i](MeteoData::HNW) = par->raincorrfact * rain
								+ par->snowcorrfact * snow;

						if (ta < 0 && hnw != IOUtils::nodata && ta != IOUtils::nodata) {
							std::cout <<std::endl<< "Test2: After correction St."<<i+1<<" HNW_with_correction:" << meteo[i](MeteoData::HNW)<<" TA:"<<ta<<"[ raincorr " << par->raincorrfact << " snowcorr " << par->snowcorrfact<< " ]"<<endl;
							//stop_execution();
						}
						/* fill the data into the temporary vector of vectors */
						vecMeteos.at(i).push_back(meteo[i]);
					}
				}

				// matteo added

//				for (size_t i = 0; i < vecMeteos.size()  ; i++) {
//					std::vector<MeteoData> meteo2=vecMeteos[i];
//					for(size_t j=0; j<meteo2.size(); j++){
//					double Iprec =meteo2[j](MeteoData::HNW);
//					double airT =meteo2[j](MeteoData::TA)- 273.15;
//					if (Iprec > 0 && airT < 0 && Iprec != IOUtils::nodata && airT != IOUtils::nodata) {
//								std::cout << "Tset3: After correction inside vecMeteos St."<<i+1<<" HNW: " << Iprec <<" TA:"<<airT<<endl;
//						}
//					}
//				}
			//	stop_execution();// matteo added

			/* Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector */
				io.add_to_cache(d1, meteo);

				/* Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector */
			//	io.push_meteo_data(IOManager::filtered, d1, d1, vecMeteos);

				// matteo: added
//				std::vector<MeteoData> meteo3;
//				io.getMeteoData(d1, meteo3);
//				for (size_t i = 0; i < meteo3.size()  ; i++) {
//							double Iprec =meteo3[i](MeteoData::HNW);
//							double airT =meteo3[i](MeteoData::TA)- 273.15;;
//							if (Iprec > 0 && airT < 0 && Iprec != IOUtils::nodata && airT != IOUtils::nodata) {
//								std::cout << "Tset4: After correction && push_meteo_data St."<<i+1<<" HNW: " << Iprec <<" TA:"<<airT<<endl;
//							}
//						}
		//	stop_execution();
				// matteo added

		/* Prepare the points */
		Coords point;
		std::vector<Coords> pointsVec; //create a vector of Coords objects

		double eastX, northY;

		/* Read the X,Y points listed in the GEOtop ListPoint file and prepare a vector to push in the MeteoIO interpolation */

		for (int i = 1; i <= par->chkpt->nrh; i++) {
			eastX = par->chkpt->co[i][ptX];
			northY = par->chkpt->co[i][ptY];
			point.setXY(eastX, northY, IOUtils::nodata);
			pointsVec.push_back(point);
		}

		std::cout << "\n[MeteoIO] Time to interpolate Point wise: "
					<< d1.toString(Date::ISO) << std::endl;

		/*Push the coordinate list in the MeteoIO interpolate method tom perform the interpolation on them */

		// interpolate Air Temperature
		io.interpolate(d1, dem, MeteoData::TA, pointsVec, resultTa);

		/* change TA values */
		for (unsigned int i = 0; i < resultTa.size(); i++) {
			if (resultTa[i] != IOUtils::nodata) {
				resultTa[i] -= 273.15; //back to celsius
			}
		}
		// interpolate RH
		io.interpolate(d1, dem, MeteoData::RH, pointsVec, resultRh);
		// interpolate Air Pressure
		io.interpolate(d1, dem, MeteoData::P, pointsVec, resultP);
		/* change P values */
		for (unsigned int i = 0; i < resultP.size(); i++) {
			if (resultP[i] != IOUtils::nodata) {
				resultP[i] /= 100.0;
			}
		}
		// interpolate Wind Speed
		io.interpolate(d1, dem, MeteoData::VW, pointsVec, resultVw);
		/* change VW values, if the measured value is less than the user given vmin value then set the vmin */
		for (unsigned int i = 0; i < resultVw.size(); i++) {
			if (resultVw[i] != IOUtils::nodata && resultVw[i] < par->Vmin) {
				resultVw[i] = par->Vmin; //set GEOtop vw min value
			}
		}
		// interpolate Wind Direction

		io.interpolate(d1, dem, MeteoData::DW, pointsVec, resultDw);
		// interpolate precipitation

		// matteo added
		if (hnw > 0 && ta < 0 && hnw != IOUtils::nodata && ta != IOUtils::nodata) {
			std::vector<MeteoData> meteo2;
			/* Read the meteo data for the given timestep */
			size_t numOfStations = io.getMeteoData(d1, meteo2);
			for(size_t i=0;i<numOfStations;i++){
				std::cout <<std::endl<< "Before 2D interp. St."<<i+1<<" HNW: " << meteo2[i](MeteoData::HNW)<<endl;
			}
		//	stop_execution();
		}// matteo added

		// spatial interpolation
		io.interpolate(d1, dem, MeteoData::HNW, pointsVec, resultHnw);

	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	std::cout << "[MeteoIO] Start copying point data to GEOtop format: "
			<< std::endl;

	/* Now copy all that MeteoIO interpolated data to the appropriate GEOtop grids */

	copyGridToMatrixPointWise(resultTa, met->Tgrid);
	copyGridToMatrixPointWise(resultRh, met->RHgrid);
	copyGridToMatrixPointWise(resultP, met->Pgrid);
	copyGridToMatrixPointWise(resultDw, met->Vdir);
	copyGridToMatrixPointWise(resultVw, met->Vgrid);
	copyGridToMatrixPointWise(resultHnw, wat->PrecTot);
}

void copyGridToMatrix(Grid2DObject& gridObject, DOUBLEMATRIX* myGrid) {

	/* This function copy MeteoIO Grid2DObject to the GEOtop structs */

	for (unsigned int ii = 0; ii < gridObject.nrows; ii++) {
		for (unsigned int jj = 0; jj < gridObject.ncols; jj++) {
			if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii)
					== IOUtils::nodata) {
				myGrid->co[ii + 1][jj + 1] = number_novalue;
			} else {
				myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj,
						gridObject.nrows - 1 - ii);
			}

			if (isnan(myGrid->co[ii + 1][jj + 1]) != 0) {
				cout << "Found nan in row " << ii << " column " << jj << endl;
				cout << "Original value coming from MeteoIO: row:"<<jj<< " col:"<< (gridObject.nrows - 1 - ii) << "  = "
						<< gridObject.grid2D(jj, gridObject.nrows - 1 - ii) << endl;
			}
		}
	}
}

void copyGridToMatrixPointWise(std::vector<double>& pointValues,
		DOUBLEMATRIX* myGrid) {

	/*
	 * This function copy MeteoIO Point vector to the GEOtop structs
	 * Notice that, co[1][1] is the first index to accessed, because in GEOtop there is an index offset
	 * therefor we have to start with [1][1]
	 */

	for (unsigned int i = 0; i < pointValues.size(); i++) {
		if (pointValues[i] != IOUtils::nodata) {
			myGrid->co[1][i + 1] = pointValues[i]; //co[1][1] is the first index accessed
		} else {
			myGrid->co[1][i + 1] = number_novalue; //co[1][1] is the first index accessed
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

extern "C" double ***meteoio_readMeteoData(long*** column,
		METEO_STATIONS *stations, double novalue, long nrOfVariables, PAR *par,
		TIMES *times) {
	long ncols = nrOfVariables; //the total number of meteo variables used in GEOtop (should stay fixed)

	//Date d1 holds the beginning of this simulation, d2 the end date of the simulation 
	//d1=times->time;
	//d2=times->time+par->Dt;
	Date d1 = par->init_date->co[i_sim];
	Date d2 = par->end_date->co[i_sim];
	//	Date d1((int)par->year0, 1, 1, 0, 0);
	//	d1 += par->JD0;
	//	d1 += times->time/86400; //times->time is in seconds, conversion to julian by devision
	//
	//	Date d2((int)par->year0, 1, 1, 0, 0);
	//	d2 += par->JD0;
	//	d2 += times->TH/24;      //the end of the simulation

	//Construction a BufferIOHandler and reading the meteo data through meteoio as configured in io.ini
	Config cfg("io.ini");
	IOManager io(cfg);

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

		io.getMeteoData(currentDate, vecM); //reading meteo data and meta datafor currentDate

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
		io.writeMeteoData(vecMeteo);
	} catch (std::exception& e) {
	} //Do nothing if not configured or error happens

	//Deal with meta data, that is allocate met->st and fill with data
	std::vector<StationData> vecMyStation;
	io.getStationData(d1, vecMyStation);

	initializeMetaData(vecMyStation, d1, novalue, par, stations);
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
			data[jj][ll][7] = novalue;
			data[jj][ll][8] = novalue;
			data[jj][ll][9] = novalue;
			data[jj][ll][10] = novalue;
			data[jj][ll][11] = novalue;
			data[jj][ll][12] = novalue;
			data[jj][ll][13] = novalue;
			//data[jj][ll][ncols] = end_vector;

			for (int gg = 0; gg < nrOfVariables; gg++) {
				if (data[jj][ll][gg] == IOUtils::nodata) {
					data[jj][ll][gg] = novalue;
				} else if (data[jj][ll][gg] != novalue) {
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
				if (data[jj][ll - 1][ff] != novalue)
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

	cout << "[I] MeteoIO NOVAL used          : " << novalue << endl;
	cout << "[I] MeteoIO #of meteo parameters: " << nrOfVariables << endl;

	//Testing access to the whole tensor
	for (unsigned int ii = 0; ii < vecStation.size(); ii++) {
		double d = 0.0;
		for (unsigned int ll = 0; ll < vecMeteo[ii].size(); ll++) {
			for (int kk = 0; kk < ncols; kk++) {
				//printf("%f ", data[ii][ll][kk]);
				d = data[ii][ll][kk];
			}
			//printf("\n");
		}
	}

	return data; //return pointer to heap allocated memory
}

void initializeMetaData(const std::vector<StationData>& vecStation,
		const Date& startDate, const double& novalue, PAR *par,
		METEO_STATIONS *stations) {
	unsigned int nrOfStations = vecStation.size();

	//Initialize the station data: set beginning of data
	int year, month, day, hour, minute;
	startDate.getDate(year, month, day, hour, minute);
	//	Date JD_start = startDate - Date((int)par->year0, 1, 1, 0, 0); // this is the actual elapsed simulation time

	//init station struct met->st
	stations->E = new_doublevector(nrOfStations); // East coordinate [m] of the meteo station
	stations->N = new_doublevector(nrOfStations); // North coordinate [m] of the meteo station
	stations->lat = new_doublevector(nrOfStations); // Latitude [rad] of the meteo station
	stations->lon = new_doublevector(nrOfStations); // Longitude [rad] of the meteo station
	stations->Z = new_doublevector(nrOfStations); // Elevation [m] of the meteo station
	stations->sky = new_doublevector(nrOfStations); // Sky-view-factor [-] of the meteo station
	stations->ST = new_doublevector(nrOfStations); // Standard time minus UTM [hours] of the meteo station
	stations->Vheight = new_doublevector(nrOfStations); // Wind velocity measurement height [m] (a.g.l.)
	stations->Theight = new_doublevector(nrOfStations); // Air temperature measurement height [m] (a.g.l.)
	//	stations->JD0=new_doublevector(nrOfStations);// Decimal Julian Day of the first data
	//	stations->Y0=new_longvector(nrOfStations);	// Year of the first data
	//	stations->Dt=new_doublevector(nrOfStations);	// Dt of sampling of the data [sec]
	//	stations->offset=new_longvector(nrOfStations);// offset column

	for (unsigned int ii = 1; ii <= nrOfStations; ii++) { //HACK
		std::cout << "[I] MeteoIO station " << ii << ":\n" << vecStation[ii - 1]
				<< std::endl;

		stations->E->co[ii] = vecStation[ii - 1].position.getEasting();
		stations->N->co[ii] = vecStation[ii - 1].position.getNorthing();
		stations->lat->co[ii] = vecStation[ii - 1].position.getLat() * PI
				/ 180.0; // from deg to [rad]
		stations->lon->co[ii] = vecStation[ii - 1].position.getLon() * PI
				/ 180.0; // from deg to [rad]
		stations->Z->co[ii] = vecStation[ii - 1].position.getAltitude();
		stations->sky->co[ii] = novalue;
		stations->ST->co[ii] = 0;
		stations->Vheight->co[ii] = 8.0;
		stations->Theight->co[ii] = 8.0;

	}
}

extern "C" void meteoio_interpolate_cloudiness(T_INIT* UV, PAR* par,
		double currentdate, DOUBLEMATRIX* tau_cloud_grid,
		DOUBLEVECTOR* tau_cloud_vec) {
	/*
	 This function performs the 2D interpolations of current cloudiness by calling the respective methods within MeteoIO
	 1) Data is copied from MeteoIO objects to GEOtop structs
	 2) interpolation takes place
	 4) gridded data copied back to GEOtop DOUBLEMATRIX inside the tau_cloud_grid
	 */

	Grid2DObject cloudwgrid;

	Date d1;
	d1.setMatlabDate(currentdate, TZ); // GEOtop use matlab offset of julian date

	std::cout << "\n[MeteoIO] Time to interpolate : cloudiness "
			<< d1.toString(Date::ISO) << std::endl;

	try {
		std::vector<std::vector<MeteoData> > vecMeteos;
		std::vector<MeteoData> meteo; // Intermediate storage for storing data sets for 1 timestep
		size_t numOfStations = io.getMeteoData(d1, meteo); // Read the meteo data for the given timestep
		vecMeteos.insert(vecMeteos.begin(), meteo.size(),
				std::vector<MeteoData>()); // Allocation for the vectors

		for (size_t i = 0; i < numOfStations; i++) {
			meteo[i](MeteoData::RSWR) = tau_cloud_vec->co[i];
			vecMeteos.at(i).push_back(meteo[i]); // fill the data into the vector of vectors
		}

		// Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector
		io.push_meteo_data(IOManager::filtered, d1, d1, vecMeteos);

		// if point_sim == 1 then point-wise simulation otherwise distributed simulation

		if (par->point_sim == 1) {

			//prepare the points
			Coords point;
			std::vector<Coords> pointsVec; //create a vector of Coords objects
			std::vector<double> resultCloud;
			double eastX, northY;

			for (int i = 1; i <= par->chkpt->nrh; i++) {
				eastX = par->chkpt->co[i][ptX];
				northY = par->chkpt->co[i][ptY];
				point.setXY(eastX, northY, IOUtils::nodata);
				pointsVec.push_back(point);
			}
            /* Interpolate point wise */
			io.interpolate(d1, dem, MeteoData::RSWR, pointsVec, resultCloud);
			/* Now copy all that data to the appropriate Array */
			copyGridToMatrixPointWise(resultCloud, tau_cloud_grid);
		} else {
            /*Interpolate 2D grid*/
			io.interpolate(d1, dem, MeteoData::RSWR, cloudwgrid);
			/* Now copy all that data to the appropriate grids */
			copyGridToMatrix(cloudwgrid, tau_cloud_grid);
		}
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}
}
#endif
