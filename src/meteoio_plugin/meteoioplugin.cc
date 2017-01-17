/*
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 - 31 December 2016
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1 is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to GEOtop Foundation and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model.
 Any feedback will be highly appreciated.
 
 */


#include "meteoioplugin.h"

using namespace std;
using namespace mio;

IOManager* io;
IOManager* io_prepocessing;
DEMObject dem;

/**
 * @brief Initialization for the MeteoIO plugin.
 * - Reading the DEM (necessary for several spatial interpolations algorithms)
 *
 * @param iomanager Reference to IOManager object that will be connected to static iomanager variable io
 */
void meteoio_init(mio::IOManager& iomanager)
{
	io = &iomanager; //pointer to the iomanager instantiated in geotop.cc
	io_prepocessing = NULL;

	string cfgfile = "io_it_extra.ini";
	if (IOUtils::fileExists(cfgfile)) {
		mio::Config cfg(cfgfile);
		io_prepocessing = new IOManager(cfg);
	}

	try {
		io->readDEM(dem);
	} catch (exception& e) {
		cerr << "[ERROR] MeteoIO: " << e.what() << endl;
	}
}

void meteoio_deallocate()
{
	if (io_prepocessing != NULL) {
		delete io_prepocessing;
	}
}

void meteoio_readDEM(GeoMatrix<double>& matrix)
{
	//copy DEM to topo struct
	matrix.resize(dem.getNy()+1, dem.getNx()+1);

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

void meteoio_readMap(const std::string& filename, GeoMatrix<double>& matrix)
{
	Grid2DObject temp;
	io->read2DGrid(temp, filename + ".asc");

	matrix.resize(temp.getNy()+1, temp.getNx()+1);
	copyGridToMatrix(temp, matrix);
}

void meteoio_read2DGrid(GeoMatrix<double>& myGrid, char* _filename)
{

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

		if (geotop::common::Variables::UV->U[1] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (geotop::common::Variables::UV->U[2] != gridObject.cellsize)
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (geotop::common::Variables::UV->U[3] != gridObject.llcorner.getNorthing())
			throw IOException("Inconsistencies between 2D Grids read", AT);
		else if (geotop::common::Variables::UV->U[4] != gridObject.llcorner.getEasting())
			throw IOException("Inconsistencies between 2D Grids read", AT);

		for (unsigned int ii = 0; ii < gridObject.getNy(); ii++) {
			for (unsigned int jj = 0; jj < gridObject.getNx(); jj++) {
				if (gridObject.grid2D(jj, gridObject.getNy() - 1 - ii)
						== IOUtils::nodata) {
				//	myGrid->co[ii + 1][jj + 1] = UV->V->co[2];
					myGrid[ii + 1][jj + 1] = geotop::common::Variables::UV->V[2];
				} else {
				//	myGrid->co[ii + 1][jj + 1] = gridObject.grid2D(jj,
				//			gridObject.nrows - 1 - ii);
					myGrid[ii + 1][jj + 1] = gridObject.grid2D(jj,
											gridObject.getNy() - 1 - ii);
				}
			}
		}
		std::cout << "Read 2D Grid from file '" << filename << "' : "
				<< gridObject.getNy() << "(rows) " << gridObject.getNx() << "(cols)"
				<< std::endl;
	} catch (std::exception& e) {
		std::cerr << "[ERROR] MeteoIO: " << e.what() << std::endl;
	}
}

void meteoio_writeEsriasciiMap(const string& filename, GeoMatrix<double>& gm)
{
	// Grid2DObject  gridObject;

	if(geotop::common::Variables::UV->U[1]!=geotop::common::Variables::UV->U[2]){
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",geotop::common::Variables::UV->U[2],geotop::common::Variables::UV->U[1]);
		t_error("Fatal error");
	}

	unsigned int ncols = gm.getCols()-1;
	unsigned int nrows = gm.getRows()-1;

	// gridObject.grid2D.resize(ncols, nrows, IOUtils::nodata);
	// gridObject.llcorner.setXY(geotop::common::Variables::UV->U[4],geotop::common::Variables::UV->U[3], 0);
	// gridObject.set(ncols, nrows, geotop::common::Variables::UV->U[1] ,gridObject.llcorner,gridObject.grid2D);

	double cellsize = geotop::common::Variables::UV->U[2];
	Coords llcorner(geotop::common::Variables::UV->U[4], geotop::common::Variables::UV->U[3]);
	Grid2DObject gridObject(ncols, nrows, cellsize, llcorner, IOUtils::nodata);

	//THE FOLLOWING CODE IS HACK:
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

void meteoio_writeEsriasciiVector(const std::string& filenam, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV)
{
	// Grid2DObject  gridObject;

	if(UV->U[1]!=UV->U[2]){
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

	// gridObject.grid2D.resize(nc, nr, IOUtils::nodata);
	// gridObject.llcorner.setXY(UV->U[4],UV->U[3], 0);
	// gridObject.set(nc,nr, UV->U[1] ,gridObject.llcorner,gridObject.grid2D);

	double cellsize = UV->U[2];
	Coords llcorner(UV->U[4], UV->U[3]);
	Grid2DObject gridObject(nc, nr, cellsize, llcorner, IOUtils::nodata);

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
 * is used to perform the partitions. PSUM (old HNW) and TA are read from the MeteoIO MeteoData and then convert TA 
 * to celsius GEOtop part_snow() method is called for partitioning and then rain correction factor and 
 * snow correction factor are added
 *
 * NOTE: THIS FUNCTIONALITY IS IMPLEMENTED IN METEOIO, the appropriate filters need to be configured, e.g.:
 * PSUM::filter3    = undercatch_wmo
 * PSUM::arg3    = cst 1.3 1.1
 * where the arguments of the undercatch_wmo filter are "cst {factor for snow} {factor for mixed precipitation}"
 *
 * @param par      Pointer to GEOtop Par object, holding the values for rain and snow correction factors
 * @param meteo    A vector of MeteoData objects (one for each station) at one point in time
 */
void hnw_correction(Par* par, std::vector<mio::MeteoData>& meteo)
{
	double rain=.0, snow=.0, hnw=.0, ta=.0;
	
	for (size_t ii = 0; ii < meteo.size(); ii++) {
		hnw = meteo[ii](MeteoData::PSUM);
		
		if (meteo[ii](MeteoData::TA) != IOUtils::nodata) { /* MeteoIO deals with temperature in Kelvin */
			ta = meteo[ii](MeteoData::TA) - 273.15;
		} else {
			ta = meteo[ii](MeteoData::TA);
		}

		/* check for non-zero precipitation, temperature must be available */
		if (hnw > 0 && hnw != IOUtils::nodata && ta != IOUtils::nodata) {
			/* perform pre-processing on PSUM (old HNW) using GEOtop part_snow method */
			part_snow(hnw, &rain, &snow, ta, par->T_rain, par->T_snow);

			/* adding rain correction factor and snow correction factor */
			meteo[ii](MeteoData::PSUM) = par->raincorrfact * rain + par->snowcorrfact * snow;
		}
	}
}

/**
 * @brief This function sees if an extra IO source (io_ini_extra.ini) is configured and
 * if so attempts to get the MeteoData for the current date from the extra source. Finally
 * it copies the values for PSUM (old HNW) and CloudFactor from the retrieved MeteoData into the already
 * present MeteoData and copies it to the IOManager cache.
 *
 * @param current Date for the current timestep
 * @param meteo MeteoData for the current timestep from the original source
 */
void merge_meteo_data(Date& current, std::vector<MeteoData>& meteo)
{
	if (io_prepocessing == NULL) return;

	std::vector<MeteoData> extra_meteo;
	io_prepocessing->getMeteoData(current, extra_meteo);

	for (std::vector<MeteoData>::iterator it = extra_meteo.begin(); it != extra_meteo.end(); ++it) {
		const string& stationID  = (*it).meta.stationID;
		//const size_t cloud_index = (*it).getParameterIndex("CloudFactor");
		const size_t cloud_index = MeteoData::TAU_CLD;

		//Look for the same stationID in meteo
		for (std::vector<MeteoData>::iterator it2 = meteo.begin(); it2 != meteo.end(); ++it2) {
			if (stationID == (*it2).meta.stationID) {
				//Perform merge
				if (cloud_index != IOUtils::npos) {
					size_t cld = (*it2).addParameter("CloudFactor");
					(*it2)(cld) = (*it)(cloud_index);
				}
				
				size_t hnw_matteo = (*it2).addParameter("HNW_MATTEO");
				(*it2)(hnw_matteo)  = (*it)(MeteoData::PSUM);
				
				continue;
			}
		}
	}
	/*
	for (std::vector<MeteoData>::iterator it2 = meteo.begin(); it2 != meteo.end(); ++it2) {
		cout << (*it2).toString() << endl;
	}
	*/
	io->add_to_points_cache(current, meteo);
}

double tDew(double T, double RH, double P)
{
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

void copyGridToMatrix(Grid2DObject& gridObject, GeoMatrix<double>& myGrid, double& thr_min, double& val_min, double& thr_max, double& val_max)
{
	for (size_t ii = 0; ii < gridObject.getNy(); ii++) {
		for (size_t jj = 0; jj < gridObject.getNx(); jj++) {
			if (gridObject.grid2D(jj, gridObject.getNy() - 1 - ii) == IOUtils::nodata) {
				myGrid[ii + 1][jj + 1] = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
			} else if (gridObject.grid2D(jj, gridObject.getNy() - 1 - ii) > thr_max) {

				myGrid[ii + 1][jj + 1] = val_max;

			} else if (gridObject.grid2D(jj, gridObject.getNy() - 1 - ii) < thr_min) {

				myGrid[ii + 1][jj + 1] = val_min;

			} else {

				myGrid[ii + 1][jj + 1] = gridObject.grid2D(jj, gridObject.getNy() - 1 - ii);
			}
			//std::cout<<"myGrid->co["<<ii<<"]["<<jj<<"]"<<myGrid->co[ii + 1][jj + 1] << std::endl;
		}
	}
}

/**
 * @brief This function copy MeteoIO Grid2DObject to a GEOtop GeoMatrix object 
 * @author noori
 */
void copyGridToMatrix(Grid2DObject& gridObject, GeoMatrix<double>& myGrid)
{
	for (size_t ii = 0; ii < gridObject.getNy(); ii++) {
		for (size_t jj = 0; jj < gridObject.getNx(); jj++) {
			if (gridObject.grid2D(jj, gridObject.getNy() - 1 - ii) == IOUtils::nodata) {
				myGrid[ii + 1][jj + 1] = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
			} else {
				myGrid[ii + 1][jj + 1] = gridObject.grid2D(jj, gridObject.getNy() - 1 - ii);
			}
			//std::cout<<"myGrid->co["<<ii<<"]["<<jj<<"]"<<myGrid->co[ii + 1][jj + 1] << std::endl;
		}
	}
}

/**
 * @brief This function copies a MeteoIO Point vector<double> to the GEOtop object GeoMatrix
 * Notice that, co[1][1] is the first index to accessed, because in GEOtop there is an index offset
 * therefor we have to start with [1][1]
 *
 * @author noori
 */
void copyGridToMatrixPointWise(const std::vector<double>& pointValues, GeoMatrix<double>& myGrid, double& thr_min, double& val_min, double& thr_max, double& val_max)
{
	for (unsigned int i = 0; i < pointValues.size(); i++) {
		if (pointValues[i] == IOUtils::nodata) {
			myGrid[1][i + 1] = geotop::input::gDoubleNoValue; //using the GEOtop nodata value
		} else if (pointValues[i] > thr_max) {

			myGrid[1][i + 1] = val_max;

		} else if (pointValues[i] < thr_min) {

			myGrid[1][i + 1] = val_min;

		} else {

			myGrid[1][i + 1] = pointValues[i];
		}
	}
}

/**
 * @brief This function copies a MeteoIO Point vector<double> to the GEOtop object GeoMatrix
 * Notice that, co[1][1] is the first index to accessed, because in GEOtop there is an index offset
 * therefor we have to start with [1][1]
 *
 * @author noori
 */
void copyGridToMatrixPointWise(const std::vector<double>& pointValues, GeoMatrix<double>& myGrid)
{
	for (unsigned int i = 0; i < pointValues.size(); i++) {
		if (pointValues[i] != IOUtils::nodata) {
			myGrid[1][i + 1] = pointValues[i]; //co[1][1] is the first index accessed
		} else {
			myGrid[1][i + 1] = geotop::input::gDoubleNoValue; //co[1][1] is the first index accessed
		}
	}
}

/**
 * @brief This function changes the values in the given grid to Celsius from Kelvin
 * @param g2d The Grid2DObject that shall be changed
 */
void convertToCelsius(Grid2DObject& g2d)
{
	g2d.grid2D -= GTConst::tk; //convert back to degrees Celsius, preserve nodata
}

/**
 * @brief This function changes the values in the given grid to mbar from Pascal
 * @param g2d The Grid2DObject that shall be changed
 */
void convertToMBar(Grid2DObject& g2d)
{
	g2d.grid2D /= 100.0;
}

/**
 * @brief This function changes the VW according to the GEOtop given min value
 *
 * @param g2d The Grid2DObject that shall be changed
 * @author noori
 */
void changeVWgrid(Grid2DObject& g2d, const double& vwMin)
{
	for (unsigned int ii = 0; ii < g2d.getNx(); ii++) {
		for (unsigned int jj = 0; jj < g2d.getNy(); jj++) {
			if ((g2d.grid2D(ii, jj) == IOUtils::nodata) || (g2d.grid2D(ii, jj) < vwMin)) {
				g2d.grid2D(ii, jj) = vwMin;
			}
		}
	}
}

/**
 * @brief This function changes a grid to given value, nodata points are preserved
 *
 * @param g2d The Grid2DObject that shall be changed
 * @author noori
 */
void changeGrid(Grid2DObject& g2d, const double& val)
{
	for (unsigned int ii = 0; ii < g2d.getNx(); ii++) {
		for (unsigned int jj = 0; jj < g2d.getNy(); jj++) {
			if (g2d.grid2D(ii, jj) != IOUtils::nodata)
				g2d.grid2D(ii, jj) = val;
		}
	}
}

/**
 * @brief Function that checks in a vector of MeteoData if a valid ISWR value is present
 *
 * @param vec_meteo vector of MeteoData
 * @param first_check true if this is the first check, a special warning message will be printed
 * @param A The AllData reference
 *
 * @return true if at least one station measured ISWR
 */
bool iswr_present(const std::vector<mio::MeteoData>& vec_meteo, const bool& first_check, AllData *A) 
{
	A->M->nstcloud   = 0; // first meteo station ID (1...n) to use for the cloudiness (ISWR)
	A->M->numstcloud = 0; // counter of meteo stations containing cloud info
	mio::Config cfg = io->getConfig();
	std::string plugin_type = cfg.get("METEO", "Input");
	if (plugin_type == "GEOTOP"){
	  if (vec_meteo.size() != (A->M->st->Z.size() - 1)) {
	    cerr << "[ERROR] Inconsistency in number of stations between GEOtop and MeteoIO. Aborting iswr_present calculation!" << endl;
	    return false;
	  }
	}

	for (size_t ii=0; ii<vec_meteo.size(); ii++) {
		if (vec_meteo[ii](MeteoData::ISWR) != IOUtils::nodata) { //ISWR is measured
			if (A->M->nstcloud == 0) A->M->nstcloud = ii+1; // records the first station that has cloud info

			A->M->numstcloud++;// increase the counter of meteo stations measuring cloud

			//	A->M->st->flag_SW_meteoST->co[ii] = 1;
			A->M->st->flag_SW_meteoST[ii+1] = 1; //HACK: this index should be checked
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

/**
 * @brief This function performs the 2D interpolations of current cloudiness by calling the respective methods within MeteoIO
 * 1) Data is copied from MeteoIO objects to GEOtop structs
 * 2) interpolation takes place
 * 3) gridded data copied back to GEOtop DOUBLEMATRIX inside the tau_cloud_grid
 */
void meteoio_interpolate_cloudiness(Par* par, const double& currentdate, GeoMatrix<double>& tau_cloud_grid, GeoVector<double>& tau_cloud_vec)
{
	Grid2DObject cloudwgrid;

	Date d1;
	d1.setMatlabDate(currentdate, geotop::common::Variables::TZ); // GEOtop use matlab offset of julian date

	std::cout << "\n[MeteoIO] Time to interpolate : cloudiness " << d1.toString(Date::ISO) << std::endl;

	try {

#ifndef WRF_PLUGIN
		// Matteo e Francesco March 2015
		// This piece of code is useful only when the cloudiness is read from Meteo Data and is necessary to push it into the MeteoIO Buffer.
		// TAU_CLD is a variable already present inside the WRF dataset (through a preprocessing of ISWR) and then inside the MeteoIO Buffer,
		// so in this case, this piece of code has to be avoided.
		std::vector<std::vector<MeteoData> > vecMeteos;
		std::vector<MeteoData> meteo; // Intermediate storage for storing data sets for 1 timestep
		int numOfStations = io->getMeteoData(d1, meteo); // Read the meteo data for the given timestep
		vecMeteos.insert(vecMeteos.begin(), meteo.size(), std::vector<MeteoData>()); // Allocation for the vectors

		for (int i = 0; i < numOfStations; i++) {
			meteo[i](MeteoData::TAU_CLD) = (tau_cloud_vec[i+1] == geotop::input::gDoubleNoValue ? IOUtils::nodata : tau_cloud_vec[i+1]);
			vecMeteos.at(i).push_back(meteo[i]); // fill the data into the vector of vectors
		}

		// Bypass the internal reading of MeteoData and to performed processing and interpolation on the data of the given vector
		// io->push_meteo_data(IOManager::filtered, d1, d1, vecMeteos);
		io->clear_cache();
		io->add_to_points_cache(d1, meteo);
#endif

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
				io->interpolate(d1, dem, MeteoData::TAU_CLD, pointsVec, resultCloud);
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
			try { /* Interpolate 2D grid */
				io->getMeteoData(d1, dem, MeteoData::TAU_CLD, cloudwgrid);
			} catch (std::exception& e) {
				changeGrid(cloudwgrid, 0.5);
			}
			/* Now copy all that data to the appropriate grids */
			copyGridToMatrix(cloudwgrid, tau_cloud_grid);
		}
		// io->write2DGrid(cloudwgrid, "cloudgrid.asc");
	} catch (std::exception& e) {
		std::cerr << "[ERROR] MeteoIO: " << e.what() << std::endl;
	}
}

/**
 * @brief This function performs the 2D interpolations by calling the respective methods within MeteoIO
 * 1) Data is copied from MeteoIO objects to GEOtop structs
 * 2) Interpolation takes place
 * 3) Gridded data copied back to GEOtop DOUBLEMATRIX
 * 4) TA, P and RH values need to be converted as well as nodata values
 *
 * @param par Pointer to the GEOtop Par object, holding geotop.inpts configuration
 * @param matlabdate Matlab julian date for which interpolation is desired
 * @param met Pointer to the GEOtop Meteo object 
 * @param wat Pointer to the GEOtop Water object
 */
void meteoio_interpolate(Par* par, double matlabdate, Meteo* met, Water* wat) {
	// We need some intermediate storage for storing the interpolated grid by MeteoIO
	Grid2DObject tagrid, rhgrid, pgrid, vwgrid, dwgrid, hnwgrid, cloudwgrid;
	Grid2DObject ilwrgrid;// grid for the interpolated ILWR
	Date current_date;
	current_date.setMatlabDate(matlabdate, geotop::common::Variables::TZ); // GEOtop use matlab offset of julian date 

	try {
		// Intermediate storage for storing data sets for 1 timestep
		std::vector<mio::MeteoData> meteo;

		// Read the meteo data for the given timestep
		io->getMeteoData(current_date, meteo);
		merge_meteo_data(current_date, meteo);

		//Bypass the internal reading of MeteoData and to performed processing and interpolation 
		//on the data of the given vector meteo
		//hnw_correction(par, meteo);
		//io->add_to_cache(current_date, meteo);

		io->getMeteoData(current_date, dem, MeteoData::TA, tagrid);
		convertToCelsius(tagrid);
		// io->write2DGrid(tagrid, MeteoGrids::TA, current_date);
		
		io->getMeteoData(current_date, dem, MeteoData::RH, rhgrid); //values as fractions from [0;1]
#ifdef ILWR_PRESENT
		io->getMeteoData(current_date, dem, MeteoData::ILWR, ilwrgrid); //TODO: to be added once WRF has ilwr
#endif
		// io->write2DGrid(ilwrgrid, MeteoGrids::ILWR, current_date);

		io->getMeteoData(current_date, dem, MeteoData::P, pgrid);
		convertToMBar(pgrid); //convert from Pascal to mbar

		try {
			io->getMeteoData(current_date, dem, MeteoData::VW, vwgrid);
			changeVWgrid(vwgrid, par->Vmin); 
		} catch (exception& e) {
			changeVWgrid(vwgrid, par->Vmin); //if something goes wrong, set to Vmin everywhere
		}

		try {
			io->getMeteoData(current_date, dem, MeteoData::DW, dwgrid);
		} catch (exception& e) {
		   changeGrid(dwgrid, 90); //if something goes wrong, set to 90 degrees everywhere
		}

		io->getMeteoData(current_date, dem, MeteoData::PSUM, hnwgrid);

#ifdef USE_DA_PLUGIN
		string cfgfile_datassim = "io_datassim.ini";
		if (IOUtils::fileExists(cfgfile_datassim)) {
			pseudo_datassim(cfgfile_datassim, current_date, MeteoData::PSUM, hnwgrid);
			pseudo_datassim_TA(cfgfile_datassim, current_date, hnwgrid, tagrid);

			// ostringstream ss;
			// ss << current_date.toString(Date::ISO) << "_TA_DA.asc";
			// io->write2DGrid(hnwgrid, MeteoGrids::PSUM, current_date);
			// io->write2DGrid(tagrid, ss.str());
		}
#endif
		//io.write2DGrid(vwgrid, "vw_change.asc");
	} catch (exception& e) {
		cerr << "[ERROR] MeteoIO: " << e.what() << endl;
	}

	cout << "[MeteoIO] Start copying Grid to GEOtop format: " << endl;

	double thr_min, val_min=0.0, thr_max=1000.0, val_max=1000.0; //thr_min val from ini
	mio::Config cfg_it = io->getConfig();
	std::string thr_str;
	try {
	cfg_it.getValue("PSUM_MIN_INTERPOLATED", "Interpolations2D", thr_str);
	thr_min = atof(thr_str.c_str());
	} catch (std::exception& e) {thr_min=0.0;} //set default value to zero	
	
	// Now copy all that data to the appropriate GEOtop grids
	copyGridToMatrix(tagrid, met->Tgrid);
	copyGridToMatrix(rhgrid, met->RHgrid);
	copyGridToMatrix(pgrid, met->Pgrid);
	copyGridToMatrix(vwgrid, met->Vgrid);
	copyGridToMatrix(dwgrid, met->Vdir);
	copyGridToMatrix(hnwgrid, wat->PrecTot, thr_min, val_min, thr_max, val_max);
#ifdef ILWR_PRESENT
	copyGridToMatrix(ilwrgrid, met->ILWRgrid);
#endif
}

/**
 * @brief This function performs point wise interpolations by calling the respective methods within MeteoIO
 * 1) Data is copied from MeteoIO objects to GEOtop structs
 * 2) Interpolation takes place
 * 3) Gridded data copied back to GEOtop DOUBLEMATRIX
 * 4) TA, P and RH values need to be converted as well as nodata values
 *
 * @param par Pointer to the GEOtop Par object, holding geotop.inpts configuration
 * @param currentdate Matlab julian date for which interpolation is desired
 * @param met Pointer to the GEOtop Meteo object 
 * @param wat Pointer to the GEOtop Water object
 */
void meteoio_interpolate_pointwise(Par* par, double currentdate, Meteo* met, Water* wat)
{
	std::vector<double> resultTa, resultRh, resultP, resultVw, resultDw, resultHnw;
	std::vector<double> resultILWR;
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

		//Prepare the points
		Coords point;
		double eastX, northY;
		std::vector<Coords> pointsVec; //create a vector of Coords objects

		//Read the X,Y points listed in the GEOtop ListPoint file and prepare a vector to push in the MeteoIO interpolation

		for (size_t i = 1; i < par->chkpt.getRows(); i++) { 		//	for (int i = 1; i <= par->chkpt->nrh; i++) {
			eastX = par->chkpt[i][ptX];  //eastX  = par->chkpt->co[i][ptX];
			northY = par->chkpt[i][ptY]; //northY = par->chkpt->co[i][ptY];
			point.setXY(eastX, northY, IOUtils::nodata);
			pointsVec.push_back(point);
		}

		std::cout << "\n[MeteoIO] Time to interpolate Point wise: " << d1.toString(Date::ISO) << std::endl;

		//Push the coordinate list in the MeteoIO interpolate method tom perform the interpolation on them
		io->interpolate(d1, dem, MeteoData::TA, pointsVec, resultTa);
		for (size_t i = 0; i < resultTa.size(); i++) { //change TA values
			if (resultTa[i] != IOUtils::nodata) resultTa[i] -= GTConst::tk; //back to celsius
		}

		io->interpolate(d1, dem, MeteoData::RH, pointsVec, resultRh);
#ifdef ILWR_PRESENT
		io->interpolate(d1, dem, MeteoData::ILWR, pointsVec, resultILWR);//TODO: to be added once WRF has ilwr
#endif
		io->interpolate(d1, dem, MeteoData::P, pointsVec, resultP);
		for (size_t i = 0; i < resultP.size(); i++) { //change P values
			if (resultP[i] != IOUtils::nodata) resultP[i] /= 100.0;
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

		try {
			io->interpolate(d1, dem, MeteoData::DW, pointsVec, resultDw);
		} catch (std::exception& e) {
			/* change DW values to a user given  value */
			for (size_t i = 0; i < resultDw.size(); i++) {
				if (resultDw[i] != IOUtils::nodata) {
					resultDw[i] = 90;
				}
			}
		}

		io->interpolate(d1, dem, MeteoData::PSUM, pointsVec, resultHnw);
#ifdef USE_DA_PLUGIN
		string cfgfile_datassim = "io_datassim.ini";
		if (IOUtils::fileExists(cfgfile_datassim)) {

			pseudo_datassim(cfgfile_datassim, d1, MeteoData::PSUM, pointsVec, resultHnw);
			pseudo_datassim_TA_pointwise(cfgfile_datassim, d1, pointsVec, resultHnw, resultTa);

			// pseudo_datassim(cfgfile_datassim, d1, MeteoData::ILWR, pointsVec, resultILWR);

		}
#endif
	} catch (std::exception& e) {
		std::cerr << "[ERROR] MeteoIO: " << e.what() << std::endl;
	}

	cout << "[MeteoIO] Start copying point data to GEOtop format: " << endl;
	/* Now copy all that MeteoIO interpolated data to the appropriate GEOtop grids */
	
	double thr_min, val_min=0.0, thr_max=1000.0, val_max=1000.0; //thr_min val from ini
	mio::Config cfg_it = io->getConfig();
	std::string thr_str;
	try {
	cfg_it.getValue("PSUM_MIN_INTERPOLATED", "Interpolations2D", thr_str);
	thr_min = atof(thr_str.c_str());
	} catch (std::exception& e) {thr_min=0.0;} //set default value to zero

	//cout << " -- PSUM_THR_MIN:" << thr_min << endl;
	
	copyGridToMatrixPointWise(resultTa, met->Tgrid);
	copyGridToMatrixPointWise(resultRh, met->RHgrid);
	copyGridToMatrixPointWise(resultP, met->Pgrid);
	copyGridToMatrixPointWise(resultDw, met->Vdir);
	copyGridToMatrixPointWise(resultVw, met->Vgrid);
	copyGridToMatrixPointWise(resultHnw, wat->PrecTot, thr_min, val_min, thr_max, val_max);
#ifdef ILWR_PRESENT
	copyGridToMatrixPointWise(resultILWR, met->ILWRgrid);//TODO: to be added once WRF has ilwr
#endif
}


