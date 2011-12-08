#ifdef USE_NETCDF_ONGOING
#include "../libraries/fluidturtle/turtle.h"
#include "../libraries/fluidturtle/t_utilities.h"
#include "../libraries/ascii/init.h"
//#include "ncgt_output.h"
#include "ncgt_turtle2netcdf.h"
#include "gt_symbols.h"

long add_2Dmap(int ncid, DOUBLEMATRIX *m, double time, const char* dimension_time, long counter, short reinitialize, double number_novalue){
	/* define the temporal counter*/
	/*!
	 *
	 * \param ncid -  (int) pointer to the netCDF archive file
	 * \param dimension_time
	 * \param couter - counter of the unlimited dimension
	 * \param reinitialize - short. If 1 m is reinizialiez and counter is upgraded
	 * \param number_novale - NULL
	 *
	 */



	//printf("\nprinting %s",m->name);
	//all->ncid=%d, all->P->output_surfenergy=%f, all->I->time=%f,fmod(all->I->time+all->P->Dt,all->P->output_surfenergy*3600.0)=%f",all->ncid,all->P->output_surfenergy,all->I->time, fmod(all->I->time+all->P->Dt,all->P->output_surfenergy*3600.0));

	ncgt_put_double_vs_time(time,dimension_time,counter, ncid,dimension_time);
	/* rotate map and put to netCDF */
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	ncgt_put_rotate180_y_doublematrix_vs_time(m,counter, ncid, dimension_time, NC_GEOTOP_XLAT, NC_GEOTOP_YLON);
	if(reinitialize==1) {
		initmatrix(0.0, m, m, number_novalue);
		counter++; // upgrede the counter
	}
	printf("\nsono qui2 %ld \n",counter);//stop_execution();
	/* 26.11.20011 to do list:
	 * put the attributes
	 *
	 *  */


	return counter;
}
#endif
