#include "meteoioplugin.h"
//#include "../geotop/meteo.h"
#include "../geotop/times.h"
#include <sstream>
using namespace std;
using namespace mio;
extern long i_sim;
const double TZ = 1.;

DEMObject dem;
Config
		cfg(
				"/Users/noori/Documents/workspace/GEOtopCPP/src/MeteoIO_plug/io_it.ini");
IOManager io(cfg);

extern "C" DOUBLEMATRIX *meteoio_readDEM(T_INIT** UVREF) {
	DOUBLEMATRIX *myDEM = NULL;

	T_INIT* UV = (T_INIT *) malloc(sizeof(T_INIT));

	try {

		//	DEMObject dem;
		IOHandler iohandler("io_it.ini");
		iohandler.readDEM(dem);

		myDEM = new_doublematrix(dem.nrows, dem.ncols);
		UV->V = new_doublevector(2);
		UV->U = new_doublevector(4);

		UV->V->co[1] = -1.0;
		UV->V->co[2] = -9999.0;

		UV->U->co[1] = dem.cellsize;
		UV->U->co[2] = dem.cellsize;
		UV->U->co[3] = dem.llcorner.getNorthing();
		UV->U->co[4] = dem.llcorner.getEasting();

		for (unsigned int ii = 0; ii < dem.nrows; ii++) {
			for (unsigned int jj = 0; jj < dem.ncols; jj++) {
				if (dem.grid2D(jj, dem.nrows - 1 - ii) == IOUtils::nodata) {
					myDEM->co[ii + 1][jj + 1] = UV->V->co[2];
				} else {
					myDEM->co[ii + 1][jj + 1] = dem.grid2D(jj, dem.nrows - 1
							- ii);
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
				<< gridObject.nrows << "(rows) " << gridObject.ncols
				<< "(cols)" << std::endl;
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	return myGrid;
}

extern "C" void meteoio_interpolate(T_INIT* UV, PAR* par, double currentdate,
		METEO* met, WATER* wat) {
	/* 
	 This function performs the 2D interpolations by calling the respective methods within MeteoIO
	 1) Data is copied from GEOtop structs to MeteoIO objects
	 2) DEM is reread
	 3) interpolation takes place
	 4) gridded data copied back to GEOtop DOUBLEMATRIX
	 5) TA and RH values need to be converted as well as nodata values
	 */
	DOUBLEMATRIX *buf;// bla bla...
	DEMObject dem;
	Grid2DObject tagrid, rhgrid, pgrid, vwgrid, dwgrid, hnwgrid;
	std::string coordsys = "", coordparam = "";

	//Config cfg("/Users/noori/Documents/workspace/GEOtopCPP/src/MeteoIO_plug/io_it.ini");
	//IOManager io(cfg);

	double novalue = UV->V->co[2];

	Date d1;
	d1.setMatlabDate(currentdate, TZ); // GEOtop use matlab offset of julian date

	std::cout << "[MeteoIO] Time to interpolate: " << std::endl;
	// std::cout << "[MeteoIO] Time to interpolate: JDbeg"<<currentdate<<" /n " << d1 << " == " << d1.toString(Date::ISO) << std::endl;

	std::vector<MeteoData> vecMeteo;
	std::vector<StationData> vecStation;

	try {
		cfg.getValue("COORDSYS", "Input", coordsys);
		cfg.getValue("COORDPARAM", "Input", coordparam, Config::nothrow);

		Coords coord(coordsys, coordparam);
		/*
		 for(unsigned int ii=1; ii<=met->st->Z->nh; ii++){//Build up MeteoData and StationData objects

		 //conversion of values required
		 //P=met->var[i-1][met->column[i-1][iPs]];
		 double p = met->var[ii-1][met->column[ii-1][iPs]];
		 if (p == novalue) p = IOUtils::nodata;

		 double rh = met->var[ii-1][met->column[ii-1][iRh]];
		 if (rh == novalue) {
		 rh = IOUtils::nodata;
		 } else {
		 rh = rh / 100.0; //MeteoIO needs values in [0;1]
		 }

		 double ta = met->var[ii-1][met->column[ii-1][iT]];
		 if (ta == novalue) {
		 ta = IOUtils::nodata;
		 } else {
		 ta = ta + 273.15; //MeteoIO deals with temperature in Kelvin;
		 }

		 double vw = met->var[ii-1][met->column[ii-1][iWs]];
		 if (vw == novalue) vw = IOUtils::nodata;

		 double vd = met->var[ii-1][met->column[ii-1][iWd]];
		 if (vd == novalue) vd = IOUtils::nodata;

		 double easting = met->st->E->co[ii];
		 double northing = met->st->N->co[ii];
		 double altitude = met->st->Z->co[ii];
		 coord.setXY(easting, northing, altitude);
		 StationData sd(coord);
		 vecStation.push_back(sd);

		 MeteoData md(currentdate);
		 md.ta = ta;
		 md.vw = vw;
		 md.dw = vd;
		 md.rh = rh;
		 md.p = p;

		 std::cout << "[MeteoIO] E: " << easting << std::endl;
		 std::cout << "[MeteoIO] N: " << northing << std::endl;
		 std::cout << "[MeteoIO] Altitude: " << altitude << std::endl;

		 std::cout << "[MeteoIO] P: " << p << std::endl;
		 std::cout << "[MeteoIO] RH: " << rh << std::endl;
		 std::cout << "[MeteoIO] TA: " << ta << std::endl;
		 std::cout << "[MeteoIO] VW: " << vw << std::endl;
		 std::cout << "[MeteoIO] VD: " << vd << std::endl;
		 std::cout << md.toString() << endl;

		 }
		 */
		//alternative: let MeteoIO do all the work
		//bufferediohandler.readMeteoData(currentdate, vecMeteo, vecStation);
		//	io.readMeteoData(d1, vecMeteo);


		//read DEM
		io.readDEM(dem);

		//interpolate
		io.interpolate(d1, dem, MeteoData::TA, tagrid);
		io.interpolate(d1, dem, MeteoData::RH, rhgrid);
		io.interpolate(d1, dem, MeteoData::P, pgrid);
		io.interpolate(d1, dem, MeteoData::VW, vwgrid);
		//	//	io.interpolate(d1, dem, MeteoData::DW, dwgrid);
		io.interpolate(d1, dem, MeteoData::HNW, hnwgrid);
		io.write2DGrid(tagrid, "ta.asc");
		//
		//		//convert values accordingly, necessary for TA and RH
		//		changeRHgrid(rhgrid);
		//		changeTAgrid(tagrid);
	} catch (std::exception& e) {
		std::cerr << "[E] MeteoIO: " << e.what() << std::endl;
	}

	std::cout << "[MeteoIO] Start copying data to GEOtop format: " << std::endl;
	if (par->point_sim == 1) {// point-wise simulation-> use DEM

		//Now copy all that data to the appropriate grids
		buf = new_doublematrix(1, par->chkpt->nrh);// matrix with 1 row and number of columns equal to the number of simulation points

		copyGridToMatrix(novalue, tagrid, buf, par->chkpt, true);
		printf("\nciao after copy TA\n");
			print_doublematrix_elements(buf,10);

		//	stop_execution();
		// extract the points

//TODO:
		free_doublematrix(buf);
	} else {// distributed simulation
		buf = new_doublematrix(tagrid.nrows, tagrid.ncols);
		buf = doubletens_to_doublemat(met->Tgrid, buf);
		print_doublematrix_elements(buf, 10);
		stop_execution();
		//Now copy all that data to the appropriate grids
		copyGridToMatrix(novalue, tagrid, buf, par->chkpt,false);
		buf = doubletens_to_doublemat(met->RHgrid, buf);
		copyGridToMatrix(novalue, rhgrid, buf, par->chkpt,false);
		buf = doubletens_to_doublemat(met->Pgrid, buf);
		copyGridToMatrix(novalue, pgrid, buf, par->chkpt,false);
		buf = doubletens_to_doublemat(met->Vgrid, buf);
		copyGridToMatrix(novalue, vwgrid, buf, par->chkpt,false);
		buf = doubletens_to_doublemat(met->Vdir, buf);
		copyGridToMatrix(novalue, dwgrid, buf, par->chkpt,false);
		buf = doubletens_to_doublemat(wat->PrecTot, buf);
		copyGridToMatrix(novalue, hnwgrid, buf, par->chkpt,false);
		free_doublematrix(buf);
	}
}

void copyGridToMatrix(const double& novalue, const Grid2DObject& gridObject,
		DOUBLEMATRIX* myGrid, DOUBLEMATRIX *coordinates, bool is_point_wise) {

	if (is_point_wise!=true) {
		for (unsigned int ii = 0; ii < gridObject.nrows; ii++) {
			for (unsigned int jj = 0; jj < gridObject.ncols; jj++) {
				if (gridObject.grid2D(jj, gridObject.nrows - 1 - ii)
						== IOUtils::nodata) {
					myGrid->co[ii + 1][jj + 1] = novalue;
				} else {
					myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj,
							gridObject.nrows - 1 - ii);
				}
				//std::cout << myGrid->co[ii+1][jj+1] << " ";						
			}
			//std::cout << std::endl;
		}
	} else {
		std::string coordin = "", coordinparam = "";
		cfg.getValue("COORDSYS", "Input", coordin);
		cfg.getValue("COORDPARAM", "Input", coordinparam, Config::nothrow);
		Coords point(coordin, coordinparam);
		point.copyProj(gridObject.llcorner);

		double eastX, northY;
		int i;

		for (i = 1; i <= coordinates->nrh; i++) {
			eastX = coordinates->co[i][ptX];
			northY = coordinates->co[i][ptY];
			point.setXY(eastX, northY, IOUtils::nodata);
			gridObject.gridify(point);
			double pointValue = gridObject(point.getGridI(), point.getGridJ());
//			std::cout <<"X : "<< eastX << " Y: "<<northY <<"Point value: " << pointValue << std::endl;
			myGrid->co[1][i] = pointValue- 273.15; //MeteoIO deals with temperature in Kelvin;

		}
		//stop_execution();
	}
}

void changeRHgrid(Grid2DObject& g2d) {
	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) *= 100.0;
		}
	}
}

void changeTAgrid(Grid2DObject& g2d) {
	for (unsigned int ii = 0; ii < g2d.ncols; ii++) {
		for (unsigned int jj = 0; jj < g2d.nrows; jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) -= 273.15; //back to celsius
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
	//BufferedIOHandler bufferediohandler(iohandler, cfg);
	//bufferediohandler.bufferAlways(false);
	//bufferediohandler.setBufferDuration(Date(1.0), Date(10.0));

	std::vector<std::vector<MeteoData> > vecMeteo; //the dimension of this vector will be nrOfStations
	std::vector<std::vector<StationData> > vecStation;//the dimension of this vector will be nrOfStations

	cout << "[I] MeteoIO: Fetching all data between " << d1.toString(Date::ISO)
			<< " and " << d2.toString(Date::ISO) << endl;

	/*
	 * In the following section the BufferedIOHandler is used to request data for descrete times
	 * MeteoIO will deal with accumulation, linear interpolation, filtering and cleaning of the data
	 * During the loop a discrete amount of time will be added to the loop variable (one hour)
	 */
	for (Date currentDate = d1; currentDate <= d2; currentDate += Date(1.0
			/ 24.0)) {
		cout << "[I] MeteoIO: Getting data for " << currentDate.toString(
				Date::ISO) << endl;
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
	//bufferediohandler.readStationData(d1, vecMyStation);
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
		while ((*column)[jj][ii] != 999999) {//the geotop way of ending a vector
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
	stations->lat = new_doublevector(nrOfStations);// Latitude [rad] of the meteo station
	stations->lon = new_doublevector(nrOfStations);// Longitude [rad] of the meteo station
	stations->Z = new_doublevector(nrOfStations); // Elevation [m] of the meteo station
	stations->sky = new_doublevector(nrOfStations);// Sky-view-factor [-] of the meteo station
	stations->ST = new_doublevector(nrOfStations); // Standard time minus UTM [hours] of the meteo station
	stations->Vheight = new_doublevector(nrOfStations); // Wind velocity measurement height [m] (a.g.l.)
	stations->Theight = new_doublevector(nrOfStations); // Air temperature measurement height [m] (a.g.l.)
	//	stations->JD0=new_doublevector(nrOfStations);// Decimal Julian Day of the first data
	//	stations->Y0=new_longvector(nrOfStations);	// Year of the first data
	//	stations->Dt=new_doublevector(nrOfStations);	// Dt of sampling of the data [sec]
	//	stations->offset=new_longvector(nrOfStations);// offset column

	for (unsigned int ii = 1; ii <= nrOfStations; ii++) { //HACK
		std::cout << "[I] MeteoIO station " << ii << ":\n"
				<< vecStation[ii - 1] << std::endl;

		stations->E->co[ii] = vecStation[ii - 1].position.getEasting();
		stations->N->co[ii] = vecStation[ii - 1].position.getNorthing();
		stations->lat->co[ii] = vecStation[ii - 1].position.getLat() * PI
				/ 180.0;// from deg to [rad]
		stations->lon->co[ii] = vecStation[ii - 1].position.getLon() * PI
				/ 180.0;// from deg to [rad]
		stations->Z->co[ii] = vecStation[ii - 1].position.getAltitude();
		stations->sky->co[ii] = novalue;
		stations->ST->co[ii] = 0;
		stations->Vheight->co[ii] = 8.0;
		stations->Theight->co[ii] = 8.0;
		//		stations->JD0->co[ii]     = JD_start.getJulianDate();
		//		stations->Y0->co[ii]      = year;
		//		stations->Dt->co[ii]      = 3600.0; //in seconds
		//		stations->offset->co[ii]  = 1;

		//cout << "JD beginning of data: " << stations->JD0->co[ii] << endl;
		//cout << "Y0 beginning of data: " << stations->Y0->co[ii] << endl;
		//cout << "Dt beginning of data: " << stations->Dt->co[ii] << endl;
	}
}

DOUBLEMATRIX * doubletens_to_doublemat(DOUBLETENSOR * input,
		DOUBLEMATRIX * output) {
	/*
	 * this functions transforms a doubletensor to a doublematrix
	 */
	//printf("\ninput->nrl=%ld,input->nrh=%ld,input->ncl=%ld,input->nch=%ld,input->ndl=%ld,input->ndh=%ld",input->nrl,input->nrh,input->ncl,input->nch,input->ndl,input->ndh);
	long r, c;
	for (r = input->nrl; r <= input->nrh; r++) {
		for (c = input->ncl; c <= input->nch; c++) {
			output->co[r][c] = input->co[input->ndl][r][c];
		}
	}
	return output;
}

