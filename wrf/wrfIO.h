/*
 * wrfIO.h
 *
 *  Created on: Aug 31, 2012
 *      Author: matteo
 */
#ifdef USE_NETCDF
#ifndef WRFIO_H_
#define WRFIO_H_

#include <stdlib.h>
#include <meteoio/MeteoIO.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include "../geotop/datastructs.h"

#define NOVALUE -9999.0
#define ERRCODE 2
#define ERROR_MESSAGE(e,n_function,n_ncfunction) {printf("Error in %s() function: %s",n_function,n_ncfunction); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
namespace mio{

class wrfIO{ //: public IOInterface {
	public:



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	int ERR(int);
	void interrogate_ncfile();
	void open_ncfile(const std::string &);
	void get_dimensions_ID_size_val(const std::string& dimension_time,const std::string& dimension_x,const std::string& dimension_y,const std::string& );
	void ncvar_to_Grid2D(const std::string& varname, Grid2DObject&, const unsigned int time_vec_index, const double layer_depth);
	int get_ForecastedMeteoGridToMask(const Grid2DObject& pred, const Grid2DObject& mask, Grid2DObject& forc, short type_spat_interp);
	void ncvar2D_to_GridArray2D(const unsigned int varid, Array2D<double>& m);
	void ncvar3D_to_GridArray2D(const unsigned int k, const unsigned int varid, Array2D<double>& m);
	void ncvar4D_to_GridArray2D(const unsigned int k, const unsigned int layer_num, const unsigned int varid, Array2D<double>& m);
	void invert_Y_Array2D(const Array2D<double>& gridObject, const unsigned int nrows,const unsigned int ncols, Array2D<double>& Geo);
	void copyGrid2DToGeoMatrix(const Grid2DObject& gridObject, GeoMatrix<double>& myGrid);
	void copyGeoMatrixToGrid2D(const GeoMatrix<double>& Geo, Grid2DObject& gridObject);
	int getIndexValueVector(const std::vector<double>& vec,const double val,const double eps);
	void close_ncfile();
	void getMinAbsValueAndIndexVector(const std::vector<double>& vec, double* min_val, unsigned int* index_min_val);
	void AddValueVector(const std::vector<double>& vec, const double val, std::vector<double>& res);
	void Downscale_UniformAssignment(const GeoMatrix<double>& mask_geo, const GeoMatrix<double>& pred_geo, GeoMatrix<double>& forc_geo, const double xll_mask, const double yll_mask, const double cellsize_mask);
	void get_coordinates_from_Array2D(const Array2D<double>& map,const double Xll, const double Yll,const unsigned int nrows, const unsigned ncols, const double cellsize, std::vector<double>& Xcoord, std::vector<double>& Ycoord);
	void read_WRF_maps(const std::string& nc_file, const std::string dimension_x,const std::string dimension_y,
			const std::string dimension_t,const std::string var_meteo, const Grid2DObject& dem, const double current_time,
			Grid2DObject& current_grid);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int ndims_in;/* id of the number of dimensions in the netCDF file*/
	int nvars_in;/* id of the number of variables in the netCDF file*/
	int ngatts_in;/* id of the number of global attributes in the netCDF file*/
	int unlimdimid_in;/* id of the unlimited dimensions in the netCDF file*/

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	int north_dimid;/* id of the north dimension */
	int east_dimid;/* id of the east dimension */
	int time_dimid;/* id of the time dimension */
	int soilD_dimid; // id of the soil depth dimension
	int north_varid;/* id of the north coordinate variable */
	int east_varid;/* id of the east coordinate variable */
	int time_varid;/* id of the time coordinate variable */
	int soilD_varid; // id of the soil depth coordinate variable;
	size_t lengthp_time;// length (number of elements) along the time dimension
	size_t lengthp_x;// length (number of elements) along the X dimension
	size_t lengthp_y;// length (number of elements) along the Y dimension
	size_t lengthp_z_soil;// length (number of layers) along the Z dimension
	std::vector<double>east_coord;//coordinates along the X dimension
	std::vector<double>time_coord;//values along the time dimension (double)
	//std::vector<std::string>time_line_str;//coordinates along the time dimension (string)
	std::vector<double>north_coord;//coordinates along the Y dimension
	std::vector<double>z_soil_coord;//coordinates along the soil depth dimension
	float nc_Xll;// lower left x coordinate
	float nc_Yll;// lower left y coordinate
	float nc_Xhr;// upper right x coordinate
	float nc_Yhr;// upper right y coordinate
	float nc_cellsize; // cellsize of the netCDF
	unsigned int nc_ncols; // number of colons of the NetCDF grid
	unsigned int nc_nrows; // number of rows of the NetCDF grid

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	int ISWR_varid;/* id of the incoming short wave radiation meteo variable */
	int TA_varid;/* id of the Air Temperature meteo variable */
	int WD_varid;/* id of the Wind Direction meteo variable */
	int WS_varid;/* id of the Wind Speed meteo variable */
	int HNW_varid;/* id of the Precipitation meteo variable */
	int RH_varid;/* id of the Relative Humidity meteo variable */

	int ncid;/* id of the netCDF file */


	private:
		void resample_grid_linear(const double& time1, const Grid2DObject& grid1,
							 const double& time2, const Grid2DObject& grid2,
							 const double& current_time, Grid2DObject& current_grid);

		void resample_grid_accumulate(const double& time1, const Grid2DObject& grid1,
							 const double& time2, const Grid2DObject& grid2,
							 const double& current_time, Grid2DObject& current_grid);

	Config cfg;
	Grid2DObject calc_grid;
	int retval;/* Error handling. */
	//std::string filename;/* name of the netCDF file that will be read. */

	int e;/* return value of the error function*/

	static const double tz_in; //netCDF time zone
};

}//end namespace


#endif /* WRFIO_H_ */
#endif
