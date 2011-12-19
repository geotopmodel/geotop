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





int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue){
   	/* define the temporal counter*/
	/*!
	 * \param m - (void *) instantaneous variable (can be doublematrix, doublevector, doubletensor)
	 * \param m0 - (void *) cumulated variable at the previous time step to be updated (can be doublematrix, doublevector, doubletensor)
	 * \param Dt: computational time step
	 * \param number_novale - NULL
	 *
	 */
	long r;// row index
	long c; // column index
	long l; // layer index
	//void* m1;// updated matrix
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // to be done
		break;
	case NC_GEOTOP_POINT_VAR: // to be done
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
		buf0=(DOUBLEMATRIX*)m0;
		buf=(DOUBLEMATRIX*)m;
		for (r=buf0->nrl; r<=buf0->nrh; r++){
			for (c=buf0->ncl; c<=buf0->nch; c++){
				if((buf0->co[r][c]!=novalue) || ((buf0->co[r][c]!=buf0->co[r][c]) && (novalue!=novalue))){
					buf0->co[r][c]+=buf->co[r][c]*Dt;
				}
			}
		}
		m0=(void*)buf0;
		m=(void*)buf;
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		buf0=(DOUBLETENSOR*)m0;
		buf=(DOUBLETENSOR*)m;
		for (r=buf0->nrl; r<=buf0->nrh; r++){
			for (c=buf0->ncl; c<=buf0->nch; c++){
				for (l=buf0->ndl; l<=buf0->ndh; l++){
					if((buf0->co[l][r][c]!=novalue) || ((buf0->co[l][r][c]!=buf0->co[l][r][c]) && (novalue!=novalue))){
						buf0->co[l][r][c]+=buf->co[l][r][c]*Dt;
					}
				}
			}
		}
		m0=(void*)buf0;
		m=(void*)buf;
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		t_error("uncorrect number of dimensions in ncgt_var_update");
		break;

	}
	return 0;
}



int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue){
   	/* define the temporal counter*/
	/*!
	 * \param m0 - (void *) cumulated variable at the previous time step to be updated (can be doublematrix, doublevector, doubletensor)
	 * \param number_novale - NULL
	 *
	 */
	long r;// row index
	long c; // column index
	long l; // layer index
	//void* m1;// updated matrix
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // to be done
		break;
	case NC_GEOTOP_POINT_VAR: // to be done
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
		buf0=(DOUBLEMATRIX*)m0;
		for (r=buf0->nrl; r<=buf0->nrh; r++){
			for (c=buf0->ncl; c<=buf0->nch; c++){
				if((buf0->co[r][c]!=novalue) || ((buf0->co[r][c]!=buf0->co[r][c]) && (novalue!=novalue))){
					buf0->co[r][c]=0;
				}
			}
		}
		m0=(void*)buf0;
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		buf0=(DOUBLETENSOR*)m0;
		for (r=buf0->nrl; r<=buf0->nrh; r++){
			for (c=buf0->ncl; c<=buf0->nch; c++){
				for (l=buf0->ndl; l<=buf0->ndh; l++){
					if((buf0->co[l][r][c]!=novalue) || ((buf0->co[l][r][c]!=buf0->co[l][r][c]) && (novalue!=novalue))){
						buf0->co[l][r][c]=0;
					}
				}
			}
		}
		m0=(void*)buf0;
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		t_error("uncorrect number of dimensions in ncgt_var_set_to_zero");
		break;

	}
	return 0;
}


void * new_output_var(void * m0, short nlimdim, double novalue){
   	/* define the temporal counter*/
	/*!
	 * \param m0 - (void *) cumulated variable at the previous time step to be updated (can be doublematrix, doublevector, doubletensor)
	 * \param number_novale - NULL
	 *
	 */
	//void* m1;// updated matrix
	void *out=NULL;
	DOUBLEMATRIX * out2d=NULL;
	DOUBLETENSOR *out3d=NULL;
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // to be done
		break;
	case NC_GEOTOP_POINT_VAR: // to be done
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
		//DOUBLEMATRIX *buf0=(DOUBLEMATRIX*)m0;
		out2d=new_doublematrix(m0->nrh, m0->nch);
		copy_doublematrix((DOUBLEMATRIX*)m0,out2d);
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		//DOUBLETENSOR * buf0=(DOUBLETENSOR*)m0;
		//out=(void *) new_doubletensor(buf0->ndh,buf0->nrh,buf0->nch);
		//(buf0, (DOUBLETENSOR *)out);
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		t_error("uncorrect number of dimensions in new_output_var");
		break;

	}
	return out;

}

long ncgt_add_output_var_time(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
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
	long counter_new;
	counter_new=ncgt_add_output_var(ncid, m, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter, reinitialize,
					NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE);
}
#endif
