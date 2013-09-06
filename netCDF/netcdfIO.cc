/*
 * netcdfIO.cc
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
#include "netcdfIO.h"

using namespace std;
using namespace mio;

//namespace mio {
/**
 * @page necdfio netCDF
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
const double NetCDFIO::tz_in = 0.; //GRIB time zone, always UTC

//IOManager* io;


void* NetCDFIO::ncgt_new_output_var(const void* m0, const short& nlimdim, const double& novalue, const std::string& suffix, const double& print_flag)
{
   	/* define the temporal counter*/
	/*!
	 * \param m0 - (void *) instantaneous variable (can be GeoMatrix, GeoVector, GeoTensor)
	 * \param number_novale - NULL
	 * \param suffix - suffix to be added to the variable_name
	 * \param print_flag - flag on printing option
	 * \description: allocate a new output variable
	 */
	//void* m1;// updated matrix
	void *out=NULL;
	GeoVector<double> *out1d=NULL;
	GeoMatrix<double> *out2d=NULL;
	GeoTensor<double> *out3d = NULL;

	if (print_flag>0){
		switch (nlimdim) {
		case NC_GEOTOP_0DIM_VAR: // TODO
			break;
		case NC_GEOTOP_POINT_VAR: // TODO
			break;
		case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
		case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
		case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:
		case NC_GEOTOP_Z_UNSTRUCT_MAP:
			out2d = new GeoMatrix<double>(*((GeoMatrix<double> *)m0));
			out2d->name = out2d->name + suffix;
			return ((void*)out2d);
			break;

		case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
			out3d = new GeoTensor<double>(*((GeoTensor<double> *)m0));
			out3d->name = out3d->name + suffix;
			//cout << "Size of this tensor: " << out3d->getDh() << " x " << out3d->getRh() << " x " << out3d->getCh() << endl;
			//printf("\n out3d->name=%s, check 31/12/2011: DA SISTEMARE\n",out3d->name);

			return ((void*)out3d);
			break;
	//	printf("\nsono qui1 a=%ld",a);//stop_execution();
		case NC_GEOTOP_UNSTRUCT_MAP:
			out1d=new GeoVector<double>(*((GeoVector<double> *)m0));
			out1d->name=out1d->name+suffix;
			return((void*)out1d);
			break;
		default:
			printf("\nincorrect number of dimensions in new_output_var!!!\n");
			exit(1);
			break;

		}
	}
	// initialization of the new variable out
	//ncgt_var_set_to_zero(out,nlimdim, novalue);
	return out;
}

int NetCDFIO::ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue){
   	/* define the temporal counter*/
	/*!
	 * \param m0 - (void *) cumulated variable at the previous time step to be updated (can be GeoMatrix, GeoVector, GeoTensor)
	 * \param number_novale - NULL
	 *
	 */
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // TODO
		break;
	case NC_GEOTOP_POINT_VAR: // TODO
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
	case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:
	case NC_GEOTOP_Z_UNSTRUCT_MAP:
	case NC_GEOTOP_Z_UNSTRUCT_MAP_IN_CONTROL_POINT:
		((GeoMatrix<double>*)m0)->reset(0, novalue);
		break;
	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
		((GeoTensor<double>*)m0)->reset(0, novalue);
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	case NC_GEOTOP_UNSTRUCT_MAP:
	case NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT:// doublevector
		((GeoVector<double>*)m0)->reset(0, novalue);
		break;
	default:
		printf("\nincorrect number of dimensions in ncgt_var_set_to_zero\n");
		exit(1);
		break;

	}
	return 0;
}

long NetCDFIO::ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step,
                                 short nlimdim, const std::string& dimension_time,const std::string& dimension_z,const std::string& dimension_x, const std::string& dimension_y,
                                 long counter, short reinitialize, short update, short rotate_y,
                                 double number_novalue, GeoMatrix<long>* rc, long **j_cont, long total_pixel)
                                 //, const long Nl, const long Nr, const long Nc){
{
	/*!
	 *
	 * \param ncid -  (int) pointer to the netCDF archive file
	 * \param m0 - (void *) cumulated variable reported from previous print time instant to be printed (can be doublematrix, doublevector, doubletensor)
	 * \param m - (void *) instantaneous variable to be printed (can be doublematrix, doublevector, doubletensor). Must be the same type of m0.
	 * \param dimension_time
	 * \param print_time_step - (double) printing time step
	 * \param computation_time_step - (double) computational time step
	 * \param dimension_z - (char *) vertical dimension
	 * \param dimension_y - (char *) dimension 1
	 * \param dimension_x - (char *) dimension 2
	 * \param rotate_y - (short) if 1 the y dimension is rotated
	 * \param nlimdim - (short) number of limited dimensions (time excluded)
	 * \param counter - counter of the unlimited dimension
	 * \param reinitialize - short. re-initializes and/or updates the cumulated variables
	 * \param update - short. If 1 and counter is updated
	 * \param number_novale - NULL
	 * \param rc - (long GeoMatrix) - matrix of the control points
	 * OUTPUT
	 * counter_new: updated counter at which the variable will be written at a successive time
	 */
	long counter_new=counter;
	if(print_time_step>0 && fmod(time,print_time_step)<1.E-5){
		if (m0==NULL) {
			// prints m
			counter_new=ncgt_add_output_var(ncid, m, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter_new,
			update, rotate_y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
			//, Nl, Nr, Nc);
		} else {
			// prints m (instantaneous)
			counter_new=ncgt_add_output_var(ncid, m0, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter,
			NC_GEOTOP_NOUPDATE_COUNTER_TIME, rotate_y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
			//, Nl, Nr, Nc);
			// prints m0 (cumulated) and updates counter
			counter_new=ncgt_add_output_var(ncid, m, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter_new,
						update, rotate_y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
						//, Nl, Nr, Nc);
			// set to zero m0 (cumulated)
			if(reinitialize==1)	ncgt_var_set_to_zero(m0, nlimdim, number_novalue);
		}
	}else if(print_time_step>0){
		// printing time not reached: updates cumulated variable
		if(reinitialize==1) ncgt_var_update(m, m0, computation_time_step,nlimdim, number_novalue);
	}
	return counter_new;
}

//int ncgt_put_double_vs_time(double v, const char *var_name, long k, int ncid, const char *dimension_t){
int NetCDFIO::ncgt_put_double_vs_time(double v, const std::string & var_name, long k, int ncid, const std::string &dimension_t){
	/*!
	 *\param v - (double) scalar variable
	 *\param var_name - name of the scalar variable
	 *\param k - (long) counter
	 *\param ncid - pointer to the netCDF archive
	 *\param dimension_t - (char *) name of the t dimension (number of times: UNLIMITED)
	 *
	 *\brief This function writes the variable contained in a vector referred to a particular time step 2D + time variable within a NetCDF file
	 *
	 * \author Emanuele Cordano, Matteo Dall'Amico
	 * \date November 2011
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option to Turn on netCDF-4 features.
	 *
	 */
	const char *function_name="ncgt_put_double_vs_time";
	int status=NC_NOERR;
	int dvar; /* pointer to the variable*/
	int ndim=1;
	int dim[ndim];
	int dimid_t; /* time dimensional id*/
	size_t start[ndim],count[ndim];
	//int dim_already_exists=1; /* this flag verifies the existence of the dimensions */

	status=nc_redef(ncid);
	if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
	/* putting the dimension */
	/* dimension_t */
	//EV unlimited dimension requires other function
	status=nc_inq_dimid(ncid,dimension_t.c_str(),&dimid_t);
	//status=nc_inq_unlimdim(ncid,&dimid_t);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_t.c_str(),NC_UNLIMITED, &dimid_t);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_dim");
	}

	   /* The dim array is used to pass the dimids of the dimensions of
			   the netCDF variables. Both of the netCDF variables we are
			   creating share the same four dimensions. In C, the
			   unlimited dimension must come first on the list of dimids. */
	dim[0] = dimid_t;

	status=nc_inq_varid(ncid,var_name.c_str(),&dvar);
	//int nc_inq_varid (int ncid, const char *name, int *varidp);
	if (status!=NC_NOERR) {
		status=nc_def_var(ncid,var_name.c_str(),NC_DOUBLE,ndim,dim,&dvar);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_var");
	}
	//else {
		//printf("Warning in %s (nc_inq_varid) variable %s (id: %d) already exists and will be overwritten \n",function_name,var_name,dvar);
	//}

	status=nc_enddef(ncid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_enddef");
	/* These settings tell netcdf to write one timestep of data. (The
	   setting of start[0] inside the loop below tells netCDF which
	   timestep to write.) */

	count[0]=1;
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=k;

	/* Write the pretend data to the file. Although netCDF supports
	 * reading and writing subsets of data, in this case we write all
	 * the data in one operation. */

	status=nc_put_vara_double(ncid,dvar,start,count,&v);
	//int nc_put_vara_double(int ncid, int varid, const size_t start[],const size_t count[], const double *dp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_put_vara_double");




	return 0;

}



int NetCDFIO::ncgt_put_doublematrix_vs_time(GeoMatrix<double>& m, long k, int ncid, const std::string& dimension_t,  const std::string& dimension_x, const std::string & dimension_y, short rotate_y){
	/*!
	 *\param m - (GeoMatrix *) variable to be written in the NetCDF
	 *\param ncid (int) - pointer to the netCDF file
	 *\param k        - (long) number of the level (0 based) at which the xy map is printed (in time)
	 *\param dimension_x - (char *) name of the t dimension (number of times: UNLIMITED)
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *
	 *\brief This function writes the variable contained in a doublematrix as a map referred to a particular time step 2D + time variable within a NetCDF file
	 *
	 * \author Emanuele Cordano
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 */
	const char *function_name="ncgt_putdoublematrix_vs_time";
	int status=NC_NOERR;
	int dimid_t,dimid_x,dimid_y; /* pointer to the dimension of the doublematrix;*/
	int dvar; /* pointer to the variable of the doublematrix */
	int ndim=3;
	int dim[ndim];
	size_t start[ndim],count[ndim];
	char hold;

	status=nc_redef(ncid);
	if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
	//int dim_already_exists=1; /* this flag verifies the existence of the dimensions */
	/* putting the dimension */
	/* dimension_t */
	//EV unlimited dimension requires other function
	status=nc_inq_dimid(ncid,dimension_t.c_str(),&dimid_t);//COS^ VA MA ne definisce una sola (correggre anche per vect e tens)
	//status=nc_inq_unlimdim(ncid,&dimid_t); //BOH??
	//status=nc_def_dim(ncid,dimension_t,NC_UNLIMITED, &dimid_t); //SOLO PROVa
	if (status!=NC_NOERR) {
		//status=nc_def_dim(ncid,dimension_y,NC_UNLIMITED, &dimid_t);//EV
		status=nc_def_dim(ncid,dimension_t.c_str(),NC_UNLIMITED, &dimid_t); //EV
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}

	/* dimension_y */
	status=nc_inq_dimid(ncid,dimension_y.c_str(),&dimid_y);
	if (status!=NC_NOERR) {
		//status=nc_def_dim(ncid,dimension_y,m->nrh, &dimid_y);
		status=nc_def_dim(ncid,dimension_y.c_str(), m.getRows()-1, &dimid_y);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	/* dimension_x */
	status=nc_inq_dimid(ncid,dimension_x.c_str(),&dimid_x);
	if (status!=NC_NOERR) {
		//status=nc_def_dim(ncid,dimension_x,m->nch, &dimid_x);
		status=nc_def_dim(ncid,dimension_x.c_str(),m.getCols()-1, &dimid_x);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	   /* The dim array is used to pass the dimids of the dimensions of
	           the netCDF variables. Both of the netCDF variables we are
	           creating share the same four dimensions. In C, the
	           unlimited dimension must come first on the list of dimids. */
	dim[0] = dimid_t;
	dim[1] = dimid_y;
	dim[2] = dimid_x;
	//        dimids[2] = lat_dimid;
	//        dimids[3] = lon_dimid;

//	dim[0]=dimid_y;
//	dim[1]=dimid_x;
	//status=nc_inq_varid(ncid,(char *)&m.name,&dvar);
	status=nc_inq_varid(ncid,m.name.c_str(),&dvar);
	if (status!=NC_NOERR) {
		//cout << endl << "m.name=" << m.name << " m.name.c_str()=" << m.name.c_str() << endl;
		//status=nc_def_var(ncid,( char *)(&m.name),NC_DOUBLE,ndim,dim,&dvar);
		status=nc_def_var(ncid,m.name.c_str(),NC_DOUBLE,ndim,dim,&dvar);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_var");
	}
	//else {
	//	printf("Warning in %s (nc_inq_varid) variable %s (id: %d) already exists and will be overwritten \n",function_name,m->name,dvar);
	//}

	status=nc_enddef(ncid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_enddef");
	/* These settings tell netcdf to write one timestep of data. (The
	   setting of start[0] inside the loop below tells netCDF which
	   timestep to write.) */

	count[0]=1;// time dimension
	//count[1]=m->nrh;// y dimension
	count[1]=m.getRows()-1;// y dimension
	//count[2]=m->nch;// x dimension
	count[2]=m.getCols()-1;// x dimension
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=k;// time dimension
	start[1]=0;//m->nrl;//EV// y dimension
	start[2]=0;//m->ncl;//EV// y dimension
	//cout << endl << start[0] << " " << start[1] << " " << start[2] << " " << count[0] << " " << count[1] << " " << count[2] << endl;cin.get(hold);
	/* Write the pretend data to the file. Although netCDF supports
	 * reading and writing subsets of data, in this case we write all
	 * the data in one operation. */
	//printf("\nstart[0]=%ld, start[1]=%ld, start[2]=%ld, count[0]=%ld, count[1]=%ld, count[2]=%ld",start[0], start[1], start[2],count[0], count[1], count[2]);
	double pres[count[1]][count[2]];
	//cout << "m.name= " << m.name << endl;
	//cout << "start[0]=" << start[0] << " start[1]=" << start[1] << " start[2]=" << start[2] << " count[0]=" << count[0] << " count[1]=" << count[1] << " count[2]=" << count[2] << endl;
	for(unsigned int r=1; r<=count[1]; r++){
		for(unsigned int c=1; c<=count[2]; c++){
			//TODO: in case this is the right place to put the "rotate" flag
			if(rotate_y==1){
				//cout << "pres[" << count[1]-r << "][" << c-1 << "]=" << m(r,c) << " corresponding to m(r,c) where r=" << r << " c=" << c << endl;
				pres[count[1]-r][c-1]=m(r,c);// rotate yes
			} else{
				//cout << "pres[" << r-1 << "][" << c-1 << "]=" << m(r,c) << " corresponding to m(r,c) where r=" << r << " c=" << c << endl;
				pres[r-1][c-1]=m(r,c);// rotate no
			}
			//cout << pres[r-1][c-1] << " ";
		}
			//cout << endl;
	}
	//cin.get(hold);

	//status=nc_put_vara_double(ncid,dvar,start,count,&(m->co[m->nrl][m->ncl]));
	//status=nc_put_vara_double(ncid,dvar,start,count,&m(1,1));
	status=nc_put_vara_double(ncid,dvar,start,count,&pres[0][0]);
	//TODO: HACK: ATTENTION, EXPLOIT WHAT DISCUSSED WITH THOMAS AS THE GEOMATRIX BEHAVES ACCORDING TO C++ VECTOR STANDARD
	//status=nc_put_vara_double(ncid,dvar,start,count,&m[1][1]);




//	unsigned int xx=0, yy=0, r, c;
//	//m.size(,xx);
//	double pres[yy][xx];
//	cout << endl << "print: " << m.name << "  GetRows=" << m.getRows() << " GetCols=" << m.getCols() << endl;
//	for (r=1; r<m.getRows(); r++) {
//		for (c=1; c<m.getCols(); c++) {
//			//pres[ii][jj]=m(ii+1,jj+1);
//			cout << m(r,c) << " ";
//		}
//		cout << endl;
//	}

	//cout<< endl<< "ciao1"<< endl;cin.get(hold);
	//status=nc_put_vara_double(ncid,dvar,start,count,&(pres[0][0]));
	//status=nc_put_var_double(ncid,dvar,&(pres[0][0]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_put_vara_double");

	return 0;

}


//int ncgt_put_doublevector_vs_time(DOUBLEVECTOR *v, long k, int ncid, const char *dimension_t,  const char *dimension_x){
int NetCDFIO::ncgt_put_doublevector_vs_time(GeoVector<double>&v, long k, int ncid, const std::string& dimension_t, const std::string &dimension_x){
	/*!
	 *\param v - (GeoVector) variable to be written in the NetCDF
	 *\param ncid (int) - pointer to the netCDF file
	 *\param k        - (long) number of the level (0 based) at which the vector(x) is printed
	 *\param dimension_x - (char *) name of the t dimension (number of times: UNLIMITED)
	 *\param dimension_x - (char *) name of the x dimension
	 *
	 *\brief This function writes the variable contained in a vector referred to a particular time step 2D + time variable within a NetCDF file
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option to Turn on netCDF-4 features.
	 *
	 */
	//const char *function_name="ncgt_put_doublevector_vs_time";
	string function_name="ncgt_put_doublevector_vs_time";
	int status=NC_NOERR;

	int dimid_t,dimid_x; /* pointer to the dimension of the doublevector;*/
	int dvar; /* pointer to the variable of the doublevector */
	int ndim=2;
	int dim[ndim];
	size_t start[ndim],count[ndim];
	//int dim_already_exists=1; /* this flag verifies the existence of the dimensions */

	status=nc_redef(ncid);
	//if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
	/* putting the dimension */
	/* dimension_t */
	//EV unlimited dimension requires other function
	status=nc_inq_dimid(ncid,dimension_t.c_str(),&dimid_t);
	//status=nc_inq_unlimdim(ncid,&dimid_t);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_t.c_str(),NC_UNLIMITED, &dimid_t);
		if (status!=NC_NOERR)
			{
		//	ERROR_MESSAGE(status,function_name,"nc_inq_dim");
			printf("Error: %s","nc_inq_dim");
			exit(1);
			}
	}

	/* dimension_x */
	status=nc_inq_dimid(ncid,dimension_x.c_str(),&dimid_x);
	if (status!=NC_NOERR) {
	//	status=nc_def_dim(ncid,dimension_x,v->nh, &dimid_x);
		status=nc_def_dim(ncid,dimension_x.c_str(),v.size(), &dimid_x);
		if (status!=NC_NOERR)
		{
		//	ERROR_MESSAGE(status,function_name,"nc_inq_dim");
			printf("Error: %s","nc_inq_dim");
			exit(1);
		}
		}

	/* The dim array is used to pass the dimids of the dimensions of the netCDF variables. Both of the netCDF variables we are
	 * creating share the same four dimensions. In C, the unlimited dimension must come first on the list of dimids. */
	dim[0] = dimid_t;
	dim[1] = dimid_x;

	status=nc_inq_varid(ncid,v.name.c_str(),&dvar);
	if (status!=NC_NOERR) {
		status=nc_def_var(ncid,v.name.c_str(),NC_DOUBLE,ndim,dim,&dvar);
		if (status!=NC_NOERR){
		//	ERROR_MESSAGE(status,function_name,"nc_def_var");
		    printf("Error: %s","nc_def_var");
			exit(1);
			}
	}

	count[0]=1;
	count[1]=v.size();//count[1]=v->nh;
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=k;
	start[1]=0;

	/* Write the pretend data to the file. Although netCDF supports
	 * reading and writing subsets of data, in this case we write all
	 * the data in one operation. */
	status=nc_enddef(ncid);
	//if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_enddef");

	status=nc_put_vara_double(ncid,dvar,start,count,&(v[1]));
	if (status!=NC_NOERR)
		{
	//	ERROR_MESSAGE(status,function_name,"nc_put_vara_double");
		printf("Error: %s","nc_put_vara_double");
				exit(1);
		}


return 0;

}


//int ncgt_put_doublevector_from_doublematrix_vs_time(DOUBLEMATRIX *dt,long k, int ncid, const char *dimension_t,  char *suffix,
//		const char *dimension_id, GeoMatrix<long>* rc)
int NetCDFIO::ncgt_put_doublevector_from_doublematrix_vs_time(const GeoMatrix<double> &dt,long k, int ncid, const std::string& dimension_t, const  std::string& suffix,
		const std::string& dimension_id, GeoMatrix<long>* rc){

	/*!
	 *\param dt - (DOUBLEMATRIX *) variable to be written in the NetCDF
	 *\param k (long) printing counter
	 *\param ncid (int) - pointer to the netCDF file
	 *\param suffix - (char *) suffix to add to the variable name
	 *\param dimension_id - (char *) name of the id dimension (number of row)
	 *\param rc - (DOUBLEMATRIX*) matrix containing the row and column at which elements are put into the necdf archive
	 *\brief This function write the variable contained in a doublematrix within a NetCDF file
	 *
	 * \author Emanuele Cordano, Matteo Dall'Amico
	 * \date january 2011
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 */

	unsigned int id,r,c;// indexes
	int status;
	GeoVector<double> M;//DOUBLEVECTOR *M=NULL;

	M.resize(rc->getRows()-1);//M=new_doublevector(rc->getRows()-1);
	M.name=dt.name+suffix;//M->name=join_strings((char *)dt->name,suffix);

	//for(id=M->nl;id<=M->nh;id++) {
	for(id=1;id<=M.size();id++) {
	//	r=rc->co[id][1];
		r=(*rc)[id][1];
	//	c=rc->co[id][2];
		c=(*rc)[id][2];
		M[id]=dt(r,c);
	}
	status=ncgt_put_doublevector_vs_time(M,k,ncid,dimension_t,dimension_id);
	//free_doublevector(M);
	return 0;
}

//long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
//		const char* dimension_y, long counter, short update, short rotate_y, double number_novalue, GeoMatrix<long>* rc, long **j, long total_pixel){
long NetCDFIO::ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const std::string& dimension_time,const std::string& dimension_z,const std::string& dimension_x,
		const std::string&  dimension_y, long counter, short update, short rotate_y, double number_novalue, GeoMatrix<long>* rc, long **j, long total_pixel)
//, const long Nl, const long Nr, const long Nc)
{

   	/* define the temporal counter*/
	/*!
	 *
	 * \param ncid -  (int) pointer to the netCDF archive file
	 * \param m - (void *) variable to be printed (can be doublematrix, doublevector, doubletensor)
	 * \param dimension_time
	 * \param dimension_z - (char *) vertical dimension
	 * \param dimension_y - (char *) dimension 1
	 * \param dimension_x - (char *) dimension 2
	 * \param rotate_y - (short) if 1 the y dimension is rotated
	 * \param nlimdim - (short) number of limited dimensions (time excluded)
	 * \param counter - counter of the unlimited dimension
	 * \param reinitialize - short. If 1 m is re-initialized (only for nlimdim=2)
	 * \param update - short. If 1 and counter is updated
	 * \param number_novale - NULL
	 * \param rc - (LONGMATRIX) matrix containing the rows and columns of the control points
	 *
	 */
	GeoTensor<double> tmp_mv;
	GeoTensor<double> tmp_bla;
	GeoVector<double> tmp_V;
	GeoMatrix<double> tmp_M, mv1;
	char hold;
	GeoTensor<double>*mv;//DOUBLETENSOR *mv;
	//GeoMatrix<double>* mv1;//DOUBLEMATRIX *mv1;
	GeoMatrix<double> *m1, *m2;//DOUBLEMATRIX *m1, *m2;
	GeoMatrix<double> *gm1;  //Alternet to *m1,*m2
	GeoTensor<double> *gm_tens;
	GeoVector<double>* V;//DOUBLEVECTOR *V;
	unsigned int r,c,l,i;
	long npoints;
	int status;
	ncgt_put_double_vs_time(time,dimension_time,counter, ncid,dimension_time);
	//char * function_name="ncgt_add_output_var";
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // TODO (e.g. discharge at the outlet)
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
		m1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//tmp_M.resize(m1->getCols(),m1->getRows(),NOVALUE);
		tmp_M.resize(m1->getRows(),m1->getCols(),NOVALUE);
		for(r=1; r<m1->getRows(); r++) {
			for(c=1; c<m1->getCols(); c++) {
				tmp_M(r,c) = (*m1)[r][c];
				}
			}

		ncgt_put_doublematrix_vs_time(tmp_M,counter, ncid, dimension_time, dimension_x, dimension_y, rotate_y);
		break;
	case NC_GEOTOP_3D_MAP:// 3D structured tensors like snow
		//cout<< endl << "NC_GEOTOP_3D_MAP";cin.get(hold);
		gm_tens = new GeoTensor<double>(*((GeoTensor<double> *)m));

		if (rotate_y==1){
			ncgt_put_doubletensor_vs_time(*gm_tens,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z, rotate_y);
		} else {
			ncgt_put_doubletensor_vs_time(*gm_tens,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z, rotate_y);
		}
		break;
	case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:// option to print 2D variables in control points
	case NC_GEOTOP_POINT_VAR:
		m1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//tmp_M.resize(m1->getCols(),m1->getRows(),NOVALUE);
		tmp_M.resize(m1->getRows(),m1->getCols(),NOVALUE);
		for(r=1; r<m1->getRows(); r++) {
			for(c=1; c<m1->getCols(); c++) {
				tmp_M(r,c) = (*m1)[r][c];
				}
			}
		//TODO: to be debugged
		status=ncgt_put_doublevector_from_doublematrix_vs_time(tmp_M,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
		break;

	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:// option to print 3D variables in control points
	case NC_GEOTOP_Z_POINT_VAR:
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		gm_tens = new GeoTensor<double>(*((GeoTensor<double> *)m));
		//printf("\ngm_tens->getDh()=%u,gm_tens->getRh()=%u,gm_tens->getCh()=%u\n",gm_tens->getDh(),gm_tens->getRh(),gm_tens->getCh());//stop_execution();
		tmp_mv.resize(gm_tens->getDh(), gm_tens->getRh(), gm_tens->getCh());
		for (l=1; l<gm_tens->getDh(); l++) {
			for(r=1; r<gm_tens->getRh(); r++) {
				for(c=1; c<gm_tens->getCh(); c++) {
					tmp_mv(l,r,c) = (*gm_tens)[l][r][c];
					}
				}
			}
		tmp_mv.name = gm_tens->name;

		//ncgt_put_doublematrix_from_doubletensor_vs_time(*gm_tens,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		ncgt_put_doublematrix_from_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc, rotate_y);
		break;
	/* rotate map and put to netCDF */
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	case NC_GEOTOP_Z_UNSTRUCT_MAP:
		//printf("\ntime=%lf, tensor netcdf\n",time);
//		out2d0=(DOUBLEMATRIX *)m0;
//		out2d=new_doublematrix(out2d0->nrh, out2d0->nch);
//		copy_doublematrix(out2d0,out2d);
//		out2d->name=join_strings((char *)out2d0->name,suffix);
		//m1=(DOUBLEMATRIX*)m;
		gm1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
//		m2=new_doublematrix(m1->nrh, m1->nch);
//		copy_doublematrix(m1,m2);
//		npoints=m1->nch;
		npoints = gm1->getCols() -1;
		//printf("m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
//		V = new_doublevector(npoints);
		tmp_V.resize(npoints+1);
//		mv=new_doubletensor(Nl,Nr,Nc);
		Nl;
		//cout << endl << "mat& ema: tmp_mv.getRh()=" << tmp_mv.getRh() << " tmp_mv.getCh()=" << tmp_mv.getCh() << " gm1->getRows()=" << gm1->getRows() << " gm1->getCols()=" << gm1->getCols() << " Nl=" << Nl << " Nr=" << Nr << " Nc=" << Nc  <<endl;
		//tmp_mv.resize(Nl+1, tmp_mv.getRh(), tmp_mv.getCh());
		tmp_mv.resize(Nl+1, Nr+1, Nc+1);

//		mv->ndl=m2->nrl;
//		for (l=m2->nrl; l<=m2->nrh; l++){
//			for(i=1; i<=total_pixel; i++){
//			//for(i=1; i<=npoints; i++){
//				V->co[i] = m2->co[l][i];
//			}
//			for(r=1; r<= Nr; r++){
//				for (c=1; c<=Nc; c++){
//		 			if (j[r][c] > 0) {
//						mv->co[l][r][c]=V->co[j[r][c]];
//		 			}else {
//		 				mv->co[l][r][c]=number_novalue;
//		 			}
//		 		}
//		 	}
//		 }
//		mv->name=m1->name;
		//cout  << "gm1->name: " << gm1->name << endl;
		for (l=1; l<=Nl; l++) {
			//cout  << "l=" << l << endl;
			for(i=1; i<=total_pixel; i++) {
				tmp_V[i] = (*gm1)[l][i];
				}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=Nc; c++){
					if (j[r][c] > 0) {
						tmp_mv(l,r,c) = tmp_V[j[r][c]];
					}else {
						tmp_mv(l,r,c) = number_novalue;
					}
					//cout << "tmp_mv("<<l<<","<<r<<","<<c<<")="<<tmp_mv(l,r,c) << " ";
					//cout << tmp_mv(l,r,c) << " ";
				}
				//cout << endl;
			}
			//cout << endl;
		}
		tmp_mv.name = gm1->name;
		//cin.get(hold);

		//printf("\nmv->name=%s mv->ndl=%ld",mv->name, mv->ndl);
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		if (rotate_y==1){
			//ncgt_put_rotate180_y_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			//TODO: HACK: add rotate
			//ncgt_put_rotate180_y_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z,rotate_y);
		} else {
			//ncgt_put_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z,rotate_y);
		}
		//free_doubletensor(mv);
		//free_doublematrix(m2);
		//free_doublevector(V);
		break;

	case NC_GEOTOP_Z_UNSTRUCT_MAP_IN_CONTROL_POINT:
//		printf("\ntime=%lf, tensor netcdf\n",time);
//		out2d0=(DOUBLEMATRIX *)m0;
//		out2d=new_doublematrix(out2d0->nrh, out2d0->nch);
//		copy_doublematrix(out2d0,out2d);
//		out2d->name=join_strings((char *)out2d0->name,suffix);

	//	m1=(DOUBLEMATRIX*)m;
		gm1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//m2=new_doublematrix(m1->nrh, m1->nch);
		//gm2 = gm1;
		//copy_doublematrix(m1,m2);
		//npoints=m1->nch;
		npoints = gm1->getCols() -1;
		//cout << endl << "g1 rows=" << gm1->getRows() << " gm1 cols=" << gm1->getCols() << " Nl=" << Nl << " Nr=" << Nr << " Nc=" << Nc << " total pixel=" << total_pixel << endl; cin.get(hold);
		//printf("> m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);cin.get(hold);
		//printf("> gm1->ncl=%ld, gm1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",1, gm1->getCols()-1, Nl, Nr, Nc, total_pixel);stop_execution();

		//tmp_V.resize(npoints+1);
		tmp_V.resize(npoints+1);
		tmp_mv.resize(Nl+1, Nr+1, Nc+1);

		for (l=1; l<gm1->getRows(); l++) {
			for(i=1; i<=total_pixel; i++) {
				tmp_V[i] = (*gm1)[l][i];
			}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=Nc; c++){
					if (j[r][c] > 0) {
						tmp_mv(l,r,c) = tmp_V[j[r][c]];
					}else {
						tmp_mv(l,r,c) = number_novalue;
					}
				}
			}
		}
		tmp_mv.name = gm1->name;
		//printf("name=%s, name1=%s",tmp_mv.name,gm1->name);stop_execution();
		/*
		V = new_doublevector(npoints);
		mv=new_doubletensor(Nl,Nr,Nc);
		mv->ndl=m2->nrl;

		for (l=m2->nrl; l<=m2->nrh; l++){
			for(i=1; i<=total_pixel; i++){
			//for(i=1; i<=npoints; i++){
				V->co[i] = m2->co[l][i];
			}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=Nc; c++){
					if (j[r][c] > 0) {
						mv->co[l][r][c]=V->co[j[r][c]];
					}else {
						mv->co[l][r][c]=number_novalue;
					}
				}
			}
		 }
		mv->name=m1->name;
		*/
		//printf("\nmv->name=%s mv->ndl=%ld",mv->name, mv->ndl);
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		ncgt_put_doublematrix_from_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc,rotate_y);

		//free_doubletensor(mv);
		//free_doublematrix(m2);
		//free_doublevector(V);
		break;

	case NC_GEOTOP_UNSTRUCT_MAP:
		printf("\ntime=%lf, Hello netcdf maps\n",time);
		m1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//tmp_M.resize(m1->getCols(),m1->getRows(),NOVALUE);
		tmp_M.resize(m1->getRows(),m1->getCols(),NOVALUE);
		for(r=1; r<m1->getRows(); r++) {
			for(c=1; c<m1->getCols(); c++) {
				tmp_M(r,c) = (*m1)[r][c];
			}
		}
		npoints=m1->getCols();//total pixels
		tmp_V.resize(npoints);
		//mv1.resize(Nc,Nr);
		mv1.resize(Nr,Nc);
//		m1=(DOUBLEMATRIX*)m;
//		m2=new_doublematrix(m1->nrh,m1->nch);
//		copy_doublematrix(m1,m2);
//		npoints=m1->nch;
//		//printf("m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
//		mv1=new_doublematrix(Nr,Nc);
//		V = new_doublevector(npoints);
		//long npoints=m1->nch;// total_pixel

		for(i=1; i<=total_pixel; i++){
		//for(i=1; i<=npoints; i++){
			//V->co[i] = m2->co[1][i];
			tmp_V[i] = tmp_M[1][i];
			}
		for(r=1; r<= Nr; r++){
			for (c=1; c<=Nc; c++){
				 if (j[r][c] > 0) {
					//mv1->co[r][c]=V->co[j[r][c]];
					mv1(r,c)=tmp_V[j[r][c]];
				 }else {
					// mv1->co[r][c]=number_novalue;
					 mv1(r,c)=number_novalue;
			 	}
			 }
		}
		mv1.name=m1->name;
		printf("\nmv->name=%s",mv1.name.c_str());
		//ncgt_put_doublevector_from_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
		if (rotate_y==1){
			//TODO: HACK: add rotate
			//ncgt_put_rotate180_y_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time, dimension_x, dimension_y);
		} else {
			ncgt_put_doublematrix_vs_time(mv1,counter, ncid, dimension_time, dimension_x, dimension_y,rotate_y);
		}
//		free_doublematrix(mv1);
//		free_doublematrix(m2);
//		free_doublevector(V);
		break;

	case NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT:
		printf("\ntime=%lf, Hello netcdf maps\n",time);
		m1 = new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//tmp_M.resize(m1->getCols(),m1->getRows(),NOVALUE);
		tmp_M.resize(m1->getRows(),m1->getCols(),NOVALUE);
		for(r=1; r<m1->getRows(); r++) {
			for(c=1; c<m1->getCols(); c++) {
				tmp_M(r,c) = (*m1)[r][c];
			}
		}
		npoints=m1->getCols();//total pixels
		tmp_V.resize(npoints);
		mv1.resize(Nc,Nr);
//		m1=(DOUBLEMATRIX*)m;
//		m2=new_doublematrix(m1->nrh,m1->nch);
//		copy_doublematrix(m1,m2);
//		npoints=m1->nch;
		//printf("m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
//		mv1=new_doublematrix(Nr,Nc);
//		V = new_doublevector(npoints);
		//long npoints=m1->nch;// total_pixel
		for(i=1; i<=total_pixel; i++){
		//for(i=1; i<=npoints; i++){
			//V->co[i] = m2->co[1][i];
			tmp_V[i] = tmp_M[1][i];
			}
		for(r=1; r<= Nr; r++){
			for (c=1; c<=Nc; c++){
				 if (j[r][c] > 0) {
					//mv1->co[r][c]=V->co[j[r][c]];
					mv1(r,c)=tmp_V[j[r][c]];
				 }else {
					 //mv1->co[r][c]=number_novalue;
					 mv1(r,c)=number_novalue;
				}
			 }
		}
		mv1.name=m1->name;
		printf("\nmv->name=%s",mv1.name.c_str());
		status=ncgt_put_doublevector_from_doublematrix_vs_time(mv1,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
//		free_doublematrix(mv1);
//		free_doublematrix(m2);
//		free_doublevector(V);
		break;

	default:
		printf("\nincorrect number of dimensions in ncgt_add_output_var\n");
		exit(1);
		break;

	}
	if(update==1){
		counter++; // upgrade the counter
	}
	/* 26.11.2011 to do list:
	 * put the attributes
	 *
	 *  */


	return counter;

}

int NetCDFIO::ncgt_put_doubletensor_vs_time(GeoTensor<double>&t, long k, int ncid, const std::string& dimension_t,  const std::string& dimension_x, const std::string& dimension_y, const std::string& dimension_z,short rotate_y)
{
	/*!
	 *\param t - (DOUBLETENSOR *) variable to be written in the NetCDF
	 *\param ncid (int) - pointer to the netCDF file
	 *\param k        - (long) number of the level (0 based) at which the xyz map is printed
	 *\param dimension_t - (char *) name of the t dimension (number of times: UNLIMITED)
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *\param dimension_z - (char *) name of the z dimension
	 *
	 *\brief This function writes the variable contained in a doubletensor as a map referred to a particular time step 3D + time variable within a NetCDF file
	 *
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.

	 *
	 */
	const char *function_name="ncgt_putgeotensor_vs_time";
	int status=NC_NOERR;

	int dimid_t,dimid_x,dimid_y,dimid_z ;/* pointer to the dimenson of the doubletensor;*/
	int dvar; /* pointer to the variable of the doublematrix */
	int ndim=4;
	int dim[ndim];
	size_t start[ndim],count[ndim];
	char hold;

	status=nc_redef(ncid);
	if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");
	/* putting the dimension */
	/* dimension_t */
	//EV unlimited dimension requires other function
	status=nc_inq_dimid(ncid,dimension_t.c_str(),&dimid_t);
	//status=nc_inq_unlimdim(ncid,&dimid_t);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_t.c_str(),NC_UNLIMITED, &dimid_t);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	/* dimension_z */
	status=nc_inq_dimid(ncid,dimension_z.c_str(),&dimid_z);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_z.c_str(),t.getDh(), &dimid_z);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	/* dimension_y */
	status=nc_inq_dimid(ncid,dimension_y.c_str(),&dimid_y);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_y.c_str(),t.getRh()-1, &dimid_y);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	/* dimension_x */
	status=nc_inq_dimid(ncid,dimension_x.c_str(),&dimid_x);
	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension_x.c_str(),t.getCh()-1, &dimid_x);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dim");
	}
	/* The dim array is used to pass the dimids of the dimensions of the netCDF variables. Both of the netCDF variables we are creating share the same four dimensions.
	* In C, the unlimited dimension must come first on the list of dimids. */
	dim[0] = dimid_t;
	dim[1] = dimid_z;
	dim[2] = dimid_y;
	dim[3] = dimid_x;

	status=nc_inq_varid(ncid,t.name.c_str(),&dvar);
	//printf("name=%s",t.name.c_str());stop_execution();
	if (status!=NC_NOERR) {
		status=nc_def_var(ncid,t.name.c_str(),NC_DOUBLE,ndim,dim,&dvar);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_var");
	}
//	else {
//		printf("Warning in %s (nc_inq_varid) variable %s (id: %d) already exists and will be overwritten \n",function_name,t.name.c_str(),dvar);
//	}


	status=nc_enddef(ncid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_enddef");

	/* These settings tell netcdf to write one timestep of data. (The setting of start[0] inside the loop below tells netCDF which timestep to write.) */
	count[0]=1;// time dimension
	count[1]=t.getDh()-1;// layer dimension
	count[2]=t.getRh()-1;// north dimension
	count[3]=t.getCh()-1;// east dimension
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=k;
	start[1]=0;
	start[2]=0;
	start[3]=0;

	//cout << "count[0]=" << count[0] << " count[1]=" << count[1] << " count[2]=" << count[2] << " count[3]=" << count[3] << " start[0]=" << start[0] <<  " start[1]=" << start[1] << " start[2]=" << start[2] << " start[3]=" << start[3] << endl;
	//cout  << "t.name: " << t.name << endl;
	double pres[count[1]][count[2]][count[3]];
	for(unsigned int l=1; l<=count[1]; l++){
		//cout  << "l=" << l << endl;
		for(unsigned int r=1; r<=count[2]; r++){
			for(unsigned int c=1; c<=count[3]; c++){
				//TODO: in case this is the right place to put the "rotate" flag
				if (rotate_y==1){
					pres[l-1][count[2]-r][c-1]=t(l,r,c);// rotate yes
				} else{
					pres[l-1][r-1][c-1]=t(l,r,c);// rotate no
				}
				//cout << pres[l-1][r-1][c-1] << " ";
			}
			//cout << endl;
		}
		//cout << endl;
	}
	//cin.get(hold);

		/* Write the pretend data to the file. Although netCDF supports reading and writing subsets of data, in this case we write all the data in one operation. */
	status=nc_put_vara_double(ncid,dvar,start,count,&pres[0][0][0]);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_put_vara_double");



	return 0;
}

int NetCDFIO::ncgt_put_doublematrix_from_doubletensor_vs_time(const GeoTensor<double>& dt,long k, int ncid, const std::string& dimension_t, const std::string& suffix,
		const std::string& dimension_id, const std::string& dimension_z, GeoMatrix<long>* rc, short rotate_y){


	/*!
	 *\param dt - (DOUBLETENSOR *) variable to be written in the NetCDF
	 *\param k (long) printing counter
	 *\param ncid (int) - pointer to the netCDF file
	 *\param suffix - (char *) suffix to add to the variable name
	 *\param dimension_id - (char *) name of the id dimension (number of row)
	 *\param dimension_z - (char *) name of the z dimension (number of d TBC)
	 *\param rc - (DOUBLEMATRIX*) matrix containing the row and column at which elements are put into the necdf archive
	 *\brief This function write the variable contained in a doubletensor within a NetCDF file
	 *
	 * \author Emanuele Cordano
	 * \date january 2011
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 */
	GeoMatrix<double>M;//DOUBLEMATRIX *M=NULL;
	unsigned int id,l,r,c;// indexes
	char hold;
//	M=new_doublematrix(dt->ndh,rc->nrh);
	M.resize(dt.getDh(),rc->getRows());//M=new_doublematrix(dt.getDh()-1, rc->getRows()-1);
	M.name=dt.name+suffix;//M->name=join_strings((char *)dt.name.c_str(), suffix);
	//cout << "dt.getDh()=" << dt.getDh() << " rc->getRows()=" << rc->getRows() << endl;cin.get(hold);
	//for(id=M->ncl;id<=M->nch;id++) {

	//cout << "M.name=" << M.name << endl;
	for(id=1;id< M.getCols();id++) {
		//for(l=M->nrl;l<=M->nrh;l++) {
		for(l=1;l< M.getRows();l++) {
		//	r=rc->co[id][1];
			r=(*rc)[id][1];
		//	c=rc->co[id][2];
			c=(*rc)[id][2];
			//cout << "r=" << r << " c=" << c << " l=" << l << " dt(l,r,c)=" << dt(l,r,c) << endl;
			//M->co[l][id]=dt->co[l][r][c];
			//M->co[l][id]=dt(l,r,c);
			M(l,id)=dt(l,r,c);
		}
		//cout << endl;
	}

	//cin.get(hold);printf("\nprint %s\n",M.name.c_str());for (r=1; r<M.getRows(); r++) {for (c=1; c<M.getCols(); c++) { printf("m[%d,%d]=%f ",r,c,M(r,c)); }printf("\n");}
	ncgt_put_doublematrix_vs_time(M,k,ncid,dimension_t,dimension_id,dimension_z, rotate_y);
	//free_doublematrix(M);
	return 0;
}

int NetCDFIO::ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue){
   	/* define the temporal counter*/
	/*!
	 * \param m - (void *) instantaneous variable (can be doublematrix, doublevector, doubletensor)
	 * \param m0 - (void *) cumulated variable at the previous time step to be updated (can be doublematrix, doublevector, doubletensor)
	 * \param Dt: computational time step
	 * \param number_novale - NULL
	 *
	 */
	unsigned int r;// row index
	unsigned int c; // column index
	unsigned int l; // layer index
	GeoMatrix<double>* out2d;//DOUBLEMATRIX *out2d=NULL;
	GeoMatrix<double>* out2d0;//DOUBLEMATRIX *out2d0=NULL;
	GeoTensor<double>* out3d;//DOUBLETENSOR *out3d=NULL;
	GeoTensor<double>* out3d0;//DOUBLETENSOR *out3d0=NULL;

	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // TODO
		break;
	case NC_GEOTOP_POINT_VAR: // TODO
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
	case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:
		//out2d0=(DOUBLEMATRIX*)m0;
		out2d0=new GeoMatrix<double>(*((GeoMatrix<double> *)m0));
		//out2d=(DOUBLEMATRIX*)m;
		out2d=new GeoMatrix<double>(*((GeoMatrix<double> *)m));
		//for (r=out2d0->nrl; r<=out2d0->nrh; r++){
		for (r=1; r<=out2d0->getRows(); r++){
			//for (c=out2d0->ncl; c<=out2d0->nch; c++){
			for (c=1; c<=out2d0->getCols(); c++){
				//if((out2d0->co[r][c]!=novalue) || ((out2d0->co[r][c]!=out2d0->co[r][c]) && (novalue!=novalue))){
				if(((*out2d0)[r][c]!=novalue)){
					//out2d0->co[r][c]+=out2d->co[r][c]*Dt;
					(*out2d0)[r][c]+=(*out2d)[r][c]*Dt;
				}
			}
		}
		m0=(void*)out2d0;
		m=(void*)out2d;
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
		//out3d0=(DOUBLETENSOR *)m0;
		out3d0=new GeoTensor<double>(*((GeoTensor<double> *)m0));
		//out3d=(DOUBLETENSOR *)m;
		out3d=new GeoTensor<double>(*((GeoTensor<double> *)m));
		//for (r=out3d0->nrl; r<=out3d0->nrh; r++){
		for (r=1; r<=out3d0->getRh(); r++){
			//for (c=out3d0->ncl; c<=out3d0->nch; c++){
			for (c=1; c<=out3d0->getCh(); c++){
				//for (l=out3d0->ndl; l<=out3d0->ndh; l++){
				for (l=1; l<=out3d0->getDh(); l++){
					//if((out3d0->co[l][r][c]!=novalue) || ((out3d0->co[l][r][c]!=out3d0->co[l][r][c]) && (novalue!=novalue))){
					if(((*out3d0)[l][r][c]!=novalue)){
						//out3d0->co[l][r][c]+=out3d->co[l][r][c]*Dt;
						(*out3d0)[l][r][c]+=(*out3d)[l][r][c]*Dt;
					}
				}
			}
		}
		m0=(void*)out3d0;
		m=(void*)out3d;
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		printf("\nincorrect number of dimensions in ncgt_var_update\n");
		exit(1);
		break;

	}
	return 0;
}

int NetCDFIO::ncgt_put_doublevector(std::vector<double> &v, int ncid, const std::string &vec_name, const std::string &dimension){
	/*!
	 *\param v - (standard vector <double> ) variable to be written in the NetCDF
	 *\param ncid (int) - pointer to the netCDF file
	 *\param dimension - (char *) name of the dimension
	 *
	 *\brief This function write the variable contained in a standard <double> vector within a NetCDF file
	 *
	 * \author Emanuele Cordano & Matteo Dall'Amico
	 * \date August 2012
	 *
	 * \return 0 if exit is ok, otherwise an error message.
	 *
	 */
	int status=NC_NOERR;

	int dimid; /* pointer to the dimension of the doublevector;*/
	int dvar; /* pointer to the variable of the doublevector */
	int ndim=1;
	int dim[1];
	const char *function_name="ncgt_put_doublevector";
	//int dim_already_exists=1; /* this flag verifies the existence of the dimensions */

	status=nc_redef(ncid);
	if (status!=NC_NOERR) NC_GEOTOP_ERROR_MESSAGE(status,function_name,"nc_redef");

	status=nc_inq_dimid(ncid,dimension.c_str(),&dimid);

	if (status!=NC_NOERR) {
		status=nc_def_dim(ncid,dimension.c_str(),v.size(), &dimid);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_dim");
	}

	dim[0]=dimid;
	status=nc_inq_varid(ncid,vec_name.c_str(),&dvar);
	if (status!=NC_NOERR) {
		status=nc_def_var(ncid,vec_name.c_str(),NC_DOUBLE,ndim,dim,&dvar);
		if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_def_var");
	}
//	else {
//		printf("Warning in %s (nc_inq_varid) variable %s (id: %d) already exists and will be overwritten \n",function_name,vec_name,dvar);
//	}

	status=nc_enddef(ncid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_enddef");

	 /* Write the pretend data to the file. Although netCDF supports
	    * reading and writing subsets of data, in this case we write all
	    * the data in one operation. */
	//status=nc_put_var_double(ncid,dvar,&(v->co[v->nl]));
	status=nc_put_var_double(ncid,dvar,&(v[0]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_put_var_double");
	return 0;

}

void NetCDFIO::get_coordinates_from_Array2D(const Array2D<double>& map,const double Xll, const double Yll,const unsigned int nrows, const unsigned ncols, const double cellsize, std::vector<double>& Xcoord, std::vector<double>& Ycoord)
/* This function gets the vector of Y (north) and X (east) coordinates of an Array2D Object
 * Author: Matteo Dall'Amico
 * Date: July 2012
 * */
{
	unsigned int ii;
	Xcoord.resize(ncols-1,NOVALUE);
	Ycoord.resize(nrows-1,NOVALUE);
	// get the coordinates X
	Xcoord[0]=Xll+cellsize*0.5;
	for(ii=1; ii<Xcoord.size(); ii++){
		Xcoord[ii]=Xcoord[0]+cellsize*ii;
	}
	// get the coordinates Y
	Ycoord[0]=Yll+cellsize*0.5;
	for(ii=1; ii<Ycoord.size(); ii++){
		Ycoord[ii]=Ycoord[0]+cellsize*ii;
	}

}


#endif
