/*
 * netcdfIO.h
 *
 *  Created on: Jul 12, 2012
 *      Author: matteo
 */
#ifdef USE_NETCDF
#ifndef NETCDFIO_H_
#define NETCDFIO_H_

//#include <meteoio/IOInterface.h>
//#include <meteoio/Config.h>

//#include <stdio.h>
#include <stdlib.h>
#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include "../geotop/datastructs.h"
#include "gt_symbols.h"
//#include <cmath>

#define TIME_LEN 10000
#define EAST_LEN 1000
#define NORTH_LEN 1000
#define NROWS 100
#define NCOLS 100
#define NOVALUE -9999.0
#define ERRCODE 2
#define ERROR_MESSAGE(e,n_function,n_ncfunction) {printf("Error in %s() function: %s",n_function,n_ncfunction); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
//namespace mio{

class NetCDFIO{ //: public IOInterface {
	public:

	void* ncgt_new_output_var(const void * m0, const short& nlimdim, const double& novalue, const std::string& suffix, const double& print_flag);
	int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);
	int ncgt_put_double_vs_time(double v, const std::string & var_name, long k, int ncid, const std::string &dimension_t);

//	long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step,
//	                                 short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,	const char* dimension_y,
//	                                 long counter, short reinitialize, short update, short rotate_y,
//	                                 double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j_cont, long total_pixel);

	long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step,
	                                 short nlimdim, const std::string& dimension_time,const std::string& dimension_z,const std::string& dimension_x, const std::string& dimension_y,
	                                 long counter, short reinitialize, short update, short rotate_y,
	                                 double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j_cont, long total_pixel);
	                                 //, const long Nl, const long Nr, const long Nc);

	int ncgt_put_doublematrix_vs_time(GeoMatrix<double>& m, long k, int ncid, const std::string& dimension_t,  const std::string& dimension_x, const std::string & dimension_y, short rotate_y);

	int ncgt_put_doubletensor_vs_time(GeoTensor<double>&t, long k, int ncid, const std::string& dimension_t,  const std::string& dimension_x, const std::string& dimension_y, const std::string& dimension_z, short rotate_y);
	int ncgt_put_doublevector_vs_time(GeoVector<double>&v, long k, int ncid, const std::string& dimension_t, const std::string &dimension_x);
	long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const std::string& dimension_time,const std::string& dimension_z,const std::string& dimension_x,
			const std::string& dimension_y, long counter, short update, short rotate_y, double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j, long total_pixel);
	//, const long Nl, const long Nr, const long Nc);

	int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

	int ncgt_put_doublematrix_from_doubletensor_vs_time(const GeoTensor<double>& dt,long k, int ncid, const std::string& dimension_t, const std::string& suffix,
			const std::string& dimension_id, const std::string& dimension_z, GeoMatrix<long>* rc, short rotate_y);
	int ncgt_put_doublevector_from_doublematrix_vs_time(const GeoMatrix<double> &dt,long k, int ncid, const std::string& dimension_t, const std::string& suffix,
			const std::string & dimension_id, GeoMatrix<long>* rc);
	int ncgt_put_doublevector(std::vector<double> &v, int ncid, const std::string &vec_name, const std::string &dimension);
	void get_coordinates_from_Array2D(const mio::Array2D<double>& map,const double Xll, const double Yll,const unsigned int nrows, const unsigned ncols, const double cellsize, std::vector<double>& Xcoord, std::vector<double>& Ycoord);

	private:
//	Config cfg;
//	mio::Grid2DObject calc_grid;
	int retval;/* Error handling. */
	//std::string filename;/* name of the netCDF file that will be read. */

	int e;/* return value of the error function*/



	static const double tz_in; //netCDF time zone

};

//}//end namespace

#endif /* NETCDFIO_H_ */
#endif
