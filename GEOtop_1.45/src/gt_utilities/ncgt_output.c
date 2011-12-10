#ifdef USE_NETCDF_ONGOING
#include "../libraries/fluidturtle/turtle.h"
#include "../libraries/fluidturtle/t_utilities.h"
#include "../libraries/ascii/init.h"
//#include "ncgt_output.h"
#include "ncgt_turtle2netcdf.h"
#include "gt_symbols.h"

long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
		const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double number_novalue){
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
	 *
	 */
	ncgt_put_double_vs_time(time,dimension_time,counter, ncid,dimension_time);
	//char * function_name="ncgt_add_output_var";
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // to be done
		break;
	case NC_GEOTOP_POINT_VAR: // to be done
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
		if (rotate_y==1){
			ncgt_put_rotate180_y_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time, dimension_x, dimension_y);
		} else {
			ncgt_put_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time, dimension_x, dimension_y);
		}
		if(reinitialize==1) {
			initmatrix(0.0, m, m, number_novalue);
		}
		break;
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
		ncgt_put_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time, dimension_z, dimension_x);// dimension_x is the ID dimension
		break;
	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		if (rotate_y==1){
			ncgt_put_rotate180_y_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
		} else {
			ncgt_put_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
		}
		break;
	/* rotate map and put to netCDF */
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		t_error("uncorrect number of dimensions in ncgt_add_output_var");
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
#endif
