/*
 * wrfIO.cc
 *
 *  Created on: Jul 12, 2012
 *      Author: matteo
 */
/***********************************************************************************/
/*  Copyright 2012 MOUNTAIN-EERING S.R.L.                                          */
/***********************************************************************************/
/* This file is part of MeteoIO.
 MeteoIO is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MeteoIO is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifdef USE_NETCDF
#include "wrfIO.h"

using namespace std;

namespace mio {
/**
 * @page wrf
 * @section netcdf_format Format
 * The Network Common Data Form, or netCDF, is an interface to a library of data access functions for storing and retrieving data in the form of arrays. An array is an n-dimensional (where
 * n is 0, 1, 2,...) rectangular structure containing items which all have the same data type (e.g. 8-bit character, 32-bit integer). A scalar is a 0-dimensional array.
 * see http://www.unidata.ucar.edu/software/netcdf/ for more information
 *
 * @section netcdf_units Units
 * The units are supposed to be:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m s^-1
 * - wind direction in deg (0=North, clockwise)
 * - precipitation in mm h^-1
 * - radiation in W m^-2
 *
 * @section netcdf_keywords Keywords
 * This plugin uses the following keywords:
 * - STATION#: input filename (in METEOPATH). As many meteofiles as needed may be specified
 * - METEOPATH: meteo files directory where to read/write the meteofiles; [Input] and [Output] sections
 * - METEOPARAM: output file format options (ASCII or BINARY that might be followed by GZIP)
 * - SPECIALPTSFILE: a path+file name to the a file containing grid coordinates of special points of interest (for special outputs)
 *
 *
 */

IOManager* io;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int wrfIO::ERR(int k)
/* Handle errors by printing an error message and exiting with a non-zero status. */
{
	printf("Error: %s\n", nc_strerror(k));
	return 2;
}


void wrfIO::interrogate_ncfile()
/* interrogate the netCDF file to get:
 * - # of attributes, # of variables, # of dimensions, # of unlimited dimensions
 * - range of the coordinate variable (e.g. east, north, time)
 */
{
	// get dimensions, variables, attributes, unlimited dimensions
	if ((retval = nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in,&unlimdimid_in))) ERR(retval);

}


void wrfIO::open_ncfile(const std::string& file)
/* opens netCDF file */
{
	if ((retval=nc_open(file.c_str(), NC_NOWRITE, &ncid)))
		ERR(retval);
}


void wrfIO::get_dimensions_ID_size_val(const std::string& dimension_time,const std::string& dimension_x,const std::string& dimension_y, const std::string& dimension_z_soil)
/*!
	 * \param varname  (const string) - name of the variable
	 * \param dimension_x (const string) - name of the dimension_x (east)
	 * \param dimension_y (const string) - name of the dimension_y (north)
	 * \param dimension_time (const string) - name of the dimension time
	 *
	 *type of variables admitted are: 1: NC_BYTE, 2: NC_CHAR, 3: NC_SHORT, 4: NC_INT, 5: NC_FLOAT, and 6: NC_DOUBLE
	 * \author Matteo Dall'Amico, based on Emanuele Cordano's idea
	 * \date July 2012
	 *
	 *
	 */
{
	string function_name="get_dimensions_ID_size_val";
	int status;
	int ndimsp;
	char hold;
	nc_type time_type;
	//get dimension ID
	//time dimension
	status=nc_inq_dimid(ncid,dimension_time.c_str(),&time_dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
	status=nc_inq_dimlen(ncid,time_dimid,&lengthp_time);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimlen");
	//printf("\ndim len 'time' is %d and the first element is: %lf\n",(int)lengthp_time, time_coord[0]);
	cout << endl << "dim_id of " << dimension_time << " is " << time_dimid << " and the length along this dimension is " << lengthp_time  << endl;

	// Soil Depth (Z) dimension
	if (dimension_z_soil != ""){
		status=nc_inq_dimid(ncid,dimension_z_soil.c_str(),&soilD_dimid);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
		status=nc_inq_dimlen(ncid,soilD_dimid,&lengthp_z_soil);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimlen");
		cout << endl << "dim_id of " << dimension_z_soil << " is " << soilD_dimid << " and the length along this dimension is " << lengthp_z_soil  << endl;
	}
	// East (X) dimension
	status=nc_inq_dimid(ncid,dimension_x.c_str(),&east_dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
	status=nc_inq_dimlen(ncid,east_dimid,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimlen");
	cout << endl << "dim_id of " << dimension_x << " is " << east_dimid << " and the length along this dimension is " << lengthp_x  << endl;

	// North (Y) dimension
	status=nc_inq_dimid(ncid,dimension_y.c_str(),&north_dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
	status=nc_inq_dimlen(ncid,north_dimid,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimlen");
	cout << endl << "dim_id of " << dimension_y << " is " << north_dimid << " and the length along this dimension is " << lengthp_y  << endl;

	// get Variable ID of coordinate variables
	//get variable time
	status=nc_inq_varid(ncid,dimension_time.c_str(),&time_varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_varid");

	status=nc_inq_vartype(ncid, time_varid,  &time_type);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_vartype");
	status=nc_inq_varndims(ncid, time_varid, &ndimsp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_varndims");
	cout << endl << "var_id of " << dimension_time << " is " << time_varid << " and is of type " << time_type << " and has dimensions: " << ndimsp<<endl;
	time_coord.resize(lengthp_time,NOVALUE);
	//time_line_str.resize(lengthp_time,"");

//	static size_t start[] = {0, 0}; /* start at first value */
//	static size_t count[] = {lengthp_time,0};
	status=nc_get_var_double(ncid,time_varid,&(time_coord[0]));
	//else if(time_type==2) status=nc_get_var_uchar(ncid,time_varid,&(time_line_str[0].c_str()));
//	unsigned char timeline[lengthp_time][27];
//	status = nc_get_vara_uchar(ncid, time_varid, start, count, timeline[0]);

	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_get_var_double");
	//printf("\ndim len 'time' is %d and the first element is: %lf\n",(int)lengthp_time, time_coord[0]);
	cout << endl << "dim length of dim_id " << time_dimid << " is " << lengthp_time << " and the first element is " << time_coord[0] << endl;

	//get variable East (X)
	status=nc_inq_varid(ncid,dimension_x.c_str(),&east_varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
	cout << endl << "var_id of " << dimension_x << " is " << east_varid << endl;
	east_coord.resize(lengthp_x,NOVALUE);
	status=nc_get_var_double(ncid,east_varid,&east_coord[0]);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_get_var_double");
	cout << endl << "dim length of dim_id " << east_dimid << " is " << lengthp_x << " and the first element is " << east_coord[0] << endl;

	//get variable North (Y)
	status=nc_inq_varid(ncid,dimension_y.c_str(),&north_varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
	cout << endl << "var_id of " << dimension_y << " is " << north_varid << endl;
	north_coord.resize(lengthp_y,NOVALUE);
	status=nc_get_var_double(ncid,north_varid,&north_coord[0]);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_get_var_double");
	cout << endl << "dim length of dim_id " << north_dimid << " is " << lengthp_y << " and the first element is " << north_coord[0] << endl;

	//get variable soil depth (Z)
	if (dimension_z_soil != ""){
		status=nc_inq_varid(ncid,dimension_z_soil.c_str(),&soilD_varid);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_dimid");
		cout << endl << "var_id of " << dimension_z_soil << " is " << soilD_varid << endl;
		z_soil_coord.resize(lengthp_z_soil,NOVALUE);
		status=nc_get_var_double(ncid,soilD_varid,&z_soil_coord[0]);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_get_var_double");
		cout << endl << "dim length of dim_id " << soilD_dimid << " is " << lengthp_z_soil << " and the first element is " << z_soil_coord[0] << endl;
	}

	if((int)lengthp_x>1 & (int)lengthp_y>1){// the east and north coordinates are more than one element
		nc_cellsize=north_coord[1]-north_coord[0];
		nc_Xll=east_coord[0]-nc_cellsize*0.5;
		nc_Yll=north_coord[0]-nc_cellsize*0.5;
		nc_ncols=east_coord.size();
		nc_nrows=north_coord.size();
		nc_Xhr=nc_Xll+nc_ncols*nc_cellsize;
		nc_Yhr=nc_Yll+nc_nrows*nc_cellsize;
		//cout << endl<< "nc_Yll=" << nc_Yll << " nc_Xll=" << nc_Xll <<  " nc_Xhr=" << nc_Xhr << " nc_Yhr=" << nc_Yhr << " nc_cellsize=" << nc_cellsize << endl;
	}else{
		ERROR_MESSAGE(9,function_name.c_str(),"East and North element dimension must be greater than 1 element");
		}


}


void wrfIO::ncvar_to_Grid2D(const std::string& varname, Grid2DObject& m, const unsigned int time_vec_index, const double layer_depth)
{
	/*!
	 * \param ncid (const char *) - pointer to the netcdf file
	 * \param varname  (const string) - name of the variable to import
	 * \param Grid2DObject - name of Grid2D object to populate with the values of the netCDF variable
	 * \param time_in (const float) - time value
	 *
	 * \author Matteo Dall'Amico
	 * \date July 2012
	 *
	 *\return a Grid2D object containing the values referred to the variable "varname".
	 */
	char hold;

	string function_name="ncvar_to_Grid2D";

	int status;
	int varid; /* pointer to the variable */
	int ndimsp; // number of dimensions of the variable
	int k=0;
	int lay=0;
	double eps_time=100;// sec
	double eps_z_depth=10;//mm

	// get variable ID
	status=nc_inq_varid(ncid,varname.c_str(),&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_varid");
	status=nc_inq_varndims(ncid, varid, &ndimsp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name.c_str(),"nc_inq_varndims");

	cout << endl << "var_id of " << varname.c_str() << " is " << varid << " and has " << ndimsp << " dimensions" << endl;

	switch (ndimsp) {
		case 3: // e.g. T(t,Y,X)
//			k=getIndexValueVector(time_coord,time_in,eps_time);
//			cout << endl << "index of time vector corresponding to time " << time_in << " is " << k << endl;
			ncvar3D_to_GridArray2D(time_vec_index,varid, m.grid2D);
		break;

		case 4: // e.g. T(t,z,Y,X)
//			k=getIndexValueVector(time_coord,time_in,eps_time);
			lay=getIndexValueVector(z_soil_coord,layer_depth,eps_z_depth);
//			cout << endl << "index of time vector corresponding to time " << time_in << " is " << k << endl;
			cout << endl << " layer number corresponding to soil depth " << layer_depth << " is " << lay << endl;
			//cin.get(hold);
			ncvar4D_to_GridArray2D(time_vec_index,lay,varid, m.grid2D);
			break;

		default:
			cout << endl << "incorrect number of dimensions in ncgt_var_update" << endl;exit(1);
		break;
	}



	m.llcorner.setXY(nc_Xll,nc_Yll,0);
	m.set(lengthp_x,lengthp_y,nc_cellsize,m.llcorner,m.grid2D);

}



int wrfIO::get_ForecastedMeteoGridToMask(const Grid2DObject& pred, const Grid2DObject& mask, Grid2DObject& forc, short type_spat_interp)
	/* This function receives in input a prediction map of a meteo variable (pred) and interpolates
	 * the value of ''pred'' according to the mask given by ''mask'' and assigns the values to the map ''forc''
	 */
{
	char hold;
	int status;
	double Yll_mask=mask.llcorner.getNorthing(); // xll of the reference grid (dem);
	double Xll_mask=mask.llcorner.getEasting();// yll of the reference grid (dem);
	double Xhr_mask=Xll_mask+mask.ncols*mask.cellsize;
	double Yhr_mask=Yll_mask+mask.nrows*mask.cellsize;

	GeoMatrix<double> pred_geo, mask_geo, forc_geo;// GeoMatrix to do the calculations
	copyGrid2DToGeoMatrix(pred,pred_geo);
	copyGrid2DToGeoMatrix(mask,mask_geo);


	if (Yll_mask == nc_Yll & Xll_mask == nc_Xll & Xhr_mask == nc_Xhr & Yhr_mask == nc_Yhr & mask.cellsize==nc_cellsize ){

		//cout << endl<< "Yll_mask=" << Yll_mask << " nc_Yll=" << nc_Yll << " Xll_mask=" << Xll_mask << " nc_Xll=" << nc_Xll << " Xhr_mask=" << Xhr_mask << " nc_Xhr=" << nc_Xhr << " Yhr_mask=" << Yhr_mask << " nc_Yhr=" << nc_Yhr << " cellsize_mask=" << cellsize_mask << " nc_cellsize=" << nc_cellsize << endl;
		// the forecast map coincides in size and cellsize with the mask reference map
		//TODO: HACK: changing the two doesn't make any change... mmhhh... not like it
		forc.grid2D.resize(mask.ncols, mask.nrows, IOUtils::nodata);
		//forc.grid2D.resize(mask.nrows, mask.ncols, IOUtils::nodata);
		if(type_spat_interp==1){
			forc=pred;
		}
		status=1;
	} else if (Yll_mask > nc_Yll & Xll_mask > nc_Xll & Xhr_mask < nc_Xhr & Yhr_mask < nc_Yhr){
		// the forecast map is bigger than the mask reference map and it is completely surrounding the mask reference map
		//cout << endl<< "Yll_mask=" << Yll_mask << " nc_Yll=" << nc_Yll << " Xll_mask=" << Xll_mask << " nc_Xll=" << nc_Xll << " Xhr_mask=" << Xhr_mask << " nc_Xhr=" << nc_Xhr << " Yhr_mask=" << Yhr_mask << " nc_Yhr=" << nc_Yhr << " cellsize_mask=" << mask.cellsize << " nc_cellsize=" << nc_cellsize << endl;
		forc.grid2D.resize(mask.ncols, mask.nrows, IOUtils::nodata);
		//forc.grid2D.resize(mask.nrows, mask.ncols, IOUtils::nodata);
		forc=mask;
		if(type_spat_interp==1){
			Downscale_UniformAssignment(mask_geo, pred_geo, forc_geo, Xll_mask, Yll_mask, mask.cellsize);
		}
		copyGeoMatrixToGrid2D(forc_geo, forc);
		status=2;


	}else{
		//cout << endl<< "Yll_mask=" << Yll_mask << " nc_Yll=" << nc_Yll << " Xll_mask=" << Xll_mask << " nc_Xll=" << nc_Xll << " Xhr_mask=" << Xhr_mask << " nc_Xhr=" << nc_Xhr << " Yhr_mask=" << Yhr_mask << " nc_Yhr=" << nc_Yhr << " cellsize_mask=" << mask.cellsize << " nc_cellsize=" << nc_cellsize << endl;
		printf("The meteo-prediction map must be bigger and surrounding the calculation grid:Error!");
		status=3;
		exit(1);
	}
	return status;
}

void wrfIO::ncvar2D_to_GridArray2D(const unsigned int varid, Array2D<double>& m)
{
	//THOMAS METHOD
	double pres[lengthp_y][lengthp_x];
	//cout << endl << "Trying to read directly into grid: " << endl;
	m.resize(lengthp_x, lengthp_y, NOVALUE);
	unsigned int xx=0, yy=0, ii, jj;
	m.size(xx, yy);
	//cout << "Grid: " << xx << " x " << yy << endl;
	static size_t start[] = {0, 0};
	static size_t count[] = {lengthp_y, lengthp_x};
	int status;
	status=nc_get_vara_double(ncid,varid,start,count,&pres[0][0]);
	for (jj=0; jj<yy; jj++) {
		for (ii=0; ii<xx; ii++) {
			m(ii,jj)=pres[jj][ii];
			//cout << m(ii,jj) << " ";
		}
		//cout << endl;
		}
}

void wrfIO::ncvar3D_to_GridArray2D(const unsigned int k, const unsigned int varid, Array2D<double>& m)
{
	//THOMAS METHOD
	double pres[lengthp_time][lengthp_y][lengthp_x];
	//cout << endl << "Trying to read directly into grid: " << endl;
	m.resize(lengthp_x, lengthp_y, NOVALUE);
	unsigned int xx=0, yy=0;
	m.size(xx, yy);
	//cout << "Grid: " << xx << " x " << yy << endl;
	static size_t start[] = {k,0, 0};
	static size_t count[] = {1,lengthp_y, lengthp_x};
	int status;
	//status=nc_get_var_double(ncid,varid,&m.grid2D(0,0));// indexes are saved wrong
	//vector<double> myvec;
	//myvec.resize(lengthp_y*lengthp_x);
	//status=nc_get_vara_double(ncid,varid,start,count,&myvec[0]);
	//status=nc_get_vara_double(ncid,varid,start,count,&m.grid2D(0,0));// original
	status=nc_get_vara_double(ncid,varid,start,count,&pres[0][0][0]);
	for (unsigned int jj=0; jj<yy; jj++) {
		for (unsigned int ii=0; ii<xx; ii++) {
			m(ii,jj)=pres[0][jj][ii];
			//cout << m(ii,jj) << " ";
		}
		//cout << endl;
		}
}


void wrfIO::ncvar4D_to_GridArray2D(const unsigned int k, const unsigned int layer_num, const unsigned int varid, Array2D<double>& m)
{
	//THOMAS METHOD
	double pres[lengthp_time][lengthp_z_soil][lengthp_y][lengthp_x];
	//cout << endl << "Trying to read directly into grid: " << endl;
	m.resize(lengthp_x, lengthp_y, NOVALUE);
	unsigned int xx=0, yy=0;
	m.size(xx, yy);
	//cout << "Grid: " << xx << " x " << yy << endl;
	static size_t start[] = {k,layer_num,0, 0};
	static size_t count[] = {1,1,lengthp_y, lengthp_x};
	int status;
	//status=nc_get_var_double(ncid,varid,&m.grid2D(0,0));// indexes are saved wrong
	//vector<double> myvec;
	//myvec.resize(lengthp_y*lengthp_x);
	//status=nc_get_vara_double(ncid,varid,start,count,&myvec[0]);
	//status=nc_get_vara_double(ncid,varid,start,count,&m.grid2D(0,0));// original
	status=nc_get_vara_double(ncid,varid,start,count,&pres[0][0][0][0]);
	for (unsigned int jj=0; jj<yy; jj++) {
		for (unsigned int ii=0; ii<xx; ii++) {
			m(ii,jj)=pres[0][0][jj][ii];
			//cout << m(ii,jj) << " ";
		}
		//cout << endl;
		}
}



void wrfIO::invert_Y_Array2D(const Array2D<double>& orig, const unsigned int nrows,const unsigned int ncols, Array2D<double>& invert)
{
	/* This function inverts an Array2D along the Y dimension */
	unsigned int ii, jj;
	invert.resize(nrows,ncols);
	for (ii = 0; ii < nrows; ii++) {
		for (jj = 0; jj < ncols; jj++) {
			if (orig(jj, nrows - 1 - ii) == IOUtils::nodata) {
				invert[ii][jj] = NOVALUE; //using the GEOtop nodata value
			} else {
				invert[ii][jj] = orig(jj, nrows - 1 - ii);
			}
			//std::cout<<"invert->co["<<ii<<"]["<<jj<<"]"<<invert->co[ii + 1][jj + 1] << std::endl;
		}
	}
}




void wrfIO::copyGrid2DToGeoMatrix(const Grid2DObject& gridObject, GeoMatrix<double>& Geo)
{
	/* This function copies a MeteoIO Grid2DObject to the GeoMatrix according to GEOtop style */
	invert_Y_Array2D(gridObject.grid2D, gridObject.nrows,gridObject.ncols, Geo);
}

void wrfIO::copyGeoMatrixToGrid2D(const GeoMatrix<double>& Geo, Grid2DObject& gridObject)
{
	/* This function copies a GeoMatrix to a MeteoIO Grid2DObject */
	gridObject.grid2D.resize(Geo.getCols(), Geo.getRows(), IOUtils::nodata);
	for (unsigned int ii = 0; ii < Geo.getRows(); ii++) {
		for (unsigned int jj = 0; jj < Geo.getCols(); jj++) {
			//cout << "ii=" << ii << " jj=" << jj << " Geo.getRows()-1-ii=" << Geo.getRows() - 1 - ii << endl;
			if (Geo(Geo.getRows() - 1 - ii,jj) == IOUtils::nodata) {
				gridObject.grid2D(jj,ii) = IOUtils::nodata; //using the GEOtop nodata value
			} else {
				gridObject.grid2D(jj,ii) = Geo(Geo.getRows() - 1 - ii,jj);
			}
			//std::cout<<"Geo->co["<<ii<<"]["<<jj<<"]"<<Geo->co[ii + 1][jj + 1] << std::endl;
		}
	}
}

int wrfIO::getIndexValueVector(const std::vector<double>& vec,const double val, const double eps)
{
	unsigned int k;
	double delta;
	for(k=0; k<vec.size(); k++){
		delta=fabs(vec[k]-val);
		if(delta<=eps) break;
	}
	return k;
}

void wrfIO::close_ncfile()
/* closes netCDF file */
{
	if ((retval = nc_close(ncid))) ERR(retval);
}


void wrfIO::getMinAbsValueAndIndexVector(const std::vector<double>& vec, double* min_val, unsigned int* index_min_val)
/* This function calculates the minimum absolute value of a vector and the index where it occurs.
 * */
{
	unsigned int k=0;
	*min_val=fabs(vec[0]);
	//cout << "beginning: *min_val=" << *min_val << endl;
	*index_min_val=0;
	for(k=0; k<vec.size(); k++){
		if(fabs(vec[k])< *min_val) {
			*min_val=fabs(vec[k]);
			*index_min_val=k;
		}
	}
	//cout << "*min_val=" << *min_val << " *index_min_val=" << *index_min_val << endl;
}

void wrfIO::AddValueVector(const std::vector<double>& vec, const double val, std::vector<double>& res)
{
	res.resize(vec.size(),NOVALUE);
	unsigned int k=0;
	for(k=0; k<vec.size(); k++){
		res[k]=vec[k]+val;
	}
}

void wrfIO::Downscale_UniformAssignment(const GeoMatrix<double>& mask_geo, const GeoMatrix<double>& pred_geo, GeoMatrix<double>& forc_geo, const double xll_mask, const double yll_mask, const double cellsize_mask)
{
	char hold;
	unsigned int ii, jj;
	double xstar;// closest x_coordinate of the Xcoord_pred coordinate to the mask point target
	double ystar;// closest y_coordinate of the Ycoord_pred coordinate to the mask point target
	unsigned int ind_xstar;// index to closest x_coordinate of the Xcoord_pred coordinate to the mask point target
	unsigned int ind_ystar;// index to closest y_coordinate of the Ycoord_pred coordinate to the mask point target

	vector<double> Xcoord_mask;
	vector<double> Ycoord_mask;

	vector<double> buf_coordY, buf_coordX;

	forc_geo=mask_geo;

	get_coordinates_from_Array2D(mask_geo, xll_mask, yll_mask, mask_geo.getRows(), mask_geo.getCols(), cellsize_mask, Xcoord_mask, Ycoord_mask);
//	cout << "Xcoord_mask: ";for (unsigned int k=0; k<Xcoord_mask.size();k++){cout << Xcoord_mask[k] << " ";};cout << endl;
//	cout << "Ycoord_mask: ";for (unsigned int k=0; k<Ycoord_mask.size();k++){cout << Ycoord_mask[k] << " ";};cout << endl;
//	cout << endl << "mask_geo.getRows()= " << mask_geo.getRows() << " mask_geo.getCols()="<< mask_geo.getCols() << endl;
//	cout << endl << "pred_geo.getRows()= " << pred_geo.getRows() << " pred_geo.getCols()="<< pred_geo.getCols() << endl;

	// find the closest pred value and assigns it to the forc matrix
	for(ii=0; ii<mask_geo.getRows(); ii++){
		//get y_star
		AddValueVector(north_coord,-Ycoord_mask[mask_geo.getRows()-1-ii],buf_coordY);
//		cout << "Ycoord_mask[row "<<ii+1<< "]=" << Ycoord_mask[mask_geo.getRows()-1-ii] << " buf_coordY: ";for (unsigned int k=0; k<buf_coordY.size();k++){cout << -buf_coordY[buf_coordY.size()-1-k] << " ";};cout << endl;
		getMinAbsValueAndIndexVector(buf_coordY,&ystar,&ind_ystar);
		ind_ystar=pred_geo.getRows()-ind_ystar+1;// starting from the first row according to GeoMatrix style
//		cout << "ii=" << ii << " ystar=" << ystar << " ind_ystar=" << ind_ystar << endl;
		for(jj=0; jj<mask_geo.getCols(); jj++){
			//get x_star
			AddValueVector(east_coord,-Xcoord_mask[jj],buf_coordX);
			getMinAbsValueAndIndexVector(buf_coordX,&xstar,&ind_xstar);
//			cout << "jj=" << jj << " xstar=" << xstar << " ind_xstar=" << ind_xstar << endl;
//			cout << "forc_geo["<<ii+1<<"]["<<jj+1<<"]=pred_geo[" << ind_ystar+1 << "]" << "[" << ind_xstar+1 << "]=" << pred_geo(ind_ystar,ind_xstar) << endl;
			//assign value
			forc_geo(ii,jj)=pred_geo(ind_ystar,ind_xstar);
		}
	}

}

void wrfIO::get_coordinates_from_Array2D(const Array2D<double>& map,const double Xll, const double Yll,const unsigned int nrows, const unsigned ncols, const double cellsize, std::vector<double>& Xcoord, std::vector<double>& Ycoord)
/* This function gets the vector of Y (north) and X (east) coordinates of an Array2D Object
 * Author: Matteo Dall'Amico
 * Date: July 2012
 * */
{
	unsigned int ii;
	Xcoord.resize(ncols,NOVALUE);
	Ycoord.resize(nrows,NOVALUE);
	// get the coordinates X
	Xcoord[0]=Xll+cellsize*0.5;
	for(ii=1; ii<ncols; ii++){
		Xcoord[ii]=Xcoord[0]+cellsize*ii;
	}
	// get the coordinates Y
	Ycoord[0]=Yll+cellsize*0.5;
	for(ii=1; ii<nrows; ii++){
		Ycoord[ii]=Ycoord[0]+cellsize*ii;
	}

}


void wrfIO::read_WRF_maps(const std::string& nc_file, const std::string dimension_x,const std::string dimension_y,
		const std::string dimension_t,const std::string var_meteo, const Grid2DObject& dem, const double current_time,
		Grid2DObject& current_grid)
{

	Grid2DObject nc_pred1;
	Grid2DObject nc_pred2;
	Grid2DObject forecast_grid1;
	Grid2DObject forecast_grid2;
	//Grid2DObject current_grid;
	unsigned int ii,jj;
	vector<double> buf;
	double time_star;
	unsigned int ind_time_star;
	double eps_time=100;// sec
	double NODATA=-9999;
	short type_spat_interp=1; // uniform interpolation
	char hold;

	//open netCDF file
	open_ncfile(nc_file.c_str());

	// get dimensions, ID, size and values
	get_dimensions_ID_size_val(dimension_t,dimension_x, dimension_y, "");

	// Get the index of time_coord that coincides with the current time
	AddValueVector(time_coord,-current_time,buf);
	for (unsigned int k=0; k<buf.size();k++){cout << -buf[k] << " ";};cout << endl;
	getMinAbsValueAndIndexVector(buf,&time_star,&ind_time_star);
	cout << endl << "the value of time_coord most similar to " << current_time << " is " << time_star+current_time << " and the corresponding index of the vector is " << ind_time_star << endl;

	//If one of the prediction grids is within tolerance (eps_time) of current_time, then no linear re-sampling is necessary
	if(fabs(time_star-current_time)<=eps_time) {
		// the current timing coincides with the prediction map time: no need to re-sample
		//read maps at t1
		ncvar_to_Grid2D(var_meteo, nc_pred1, ind_time_star,NODATA);
		//create prediction map
		get_ForecastedMeteoGridToMask(nc_pred1, dem, forecast_grid1,type_spat_interp);
		current_grid.grid2D.resize(dem.ncols, dem.nrows, NODATA);
		current_grid=forecast_grid1;

	}else{
		// resize current_grid
		current_grid.grid2D.resize(dem.ncols, dem.nrows, NODATA);
		current_grid=dem;

		//read maps at t1 and t2
		ncvar_to_Grid2D(var_meteo, nc_pred1, ind_time_star,NODATA);
		ncvar_to_Grid2D(var_meteo, nc_pred2, ind_time_star+1,NODATA);
		//create prediction map
		get_ForecastedMeteoGridToMask(nc_pred1, dem, forecast_grid1,type_spat_interp);
		get_ForecastedMeteoGridToMask(nc_pred2, dem, forecast_grid2,type_spat_interp);

		// TODO:Check which parameter needs to be re-sampled, anything but precipitation will be linearly re-sampled
		// occorre definire il nome accettabile delle variabili (es. AirTemp, Precip, WindDir ecc.)

		if (var_meteo != "RAIN") {
		// re-sample linear
			resample_grid_linear(time_star, forecast_grid1, time_coord[ind_time_star+1], forecast_grid2, current_time, current_grid);
		}else{
		// otherwise re-sample accumulate (// aspetta funzione di Thomas)
			current_grid.grid2D.resize(current_grid.ncols,current_grid.nrows,0.0);
			//resample_grid_accumulate(time1, nc_pred1, time2, nc_pred2, current_time, current_grid);
		}

	}

//	cout << "current grid" << endl;
//	for(jj=0; jj<current_grid.ncols; jj++){
//		for(ii=0; ii<current_grid.nrows; ii++){
//			cout << current_grid(jj,ii) << " ";
//		//printf("dem[%u][%u]=%lf ",ii,jj,matr(ii,jj));
//		}
//		cout << endl;
//	}
//	cin.get(hold);

	//close the netCDF file
	close_ncfile();

}

void wrfIO::resample_grid_linear(const double& time1, const Grid2DObject& grid1,
							 const double& time2, const Grid2DObject& grid2,
							 const double& current_time, Grid2DObject& current_grid)
{
	double weight = (current_time - time1) / (time2 - time1);

	for (unsigned int jj=0; jj<grid1.nrows; jj++) {
		for (unsigned int ii=0; ii<grid1.ncols; ii++) {
			current_grid.grid2D(ii,jj) = Interpol1D::weightedMean(grid1.grid2D(ii,jj), grid2.grid2D(ii,jj), weight);
			//cout << "grid1:" << grid1.grid2D(ii,jj) << "  grid2:" << grid2.grid2D(ii,jj) << " weight:"<< weight << " current_grid:" << current_grid.grid2D(ii,jj);
		}
		//cout << endl;
	}
}

void wrfIO::resample_grid_accumulate(const double& time1, const Grid2DObject& grid1,
								const double& time2, const Grid2DObject& grid2,
								const double& current_time, Grid2DObject& current_grid)
{
	//double weight = (current_time - time1) / (time2 - time1);

	//For now just assume 0
	for (unsigned int jj=0; jj<grid1.nrows; jj++) {
		for (unsigned int ii=0; ii<grid1.ncols; ii++) {
			current_grid.grid2D(ii,jj) = 0;
		}
	}
}

}

#endif
