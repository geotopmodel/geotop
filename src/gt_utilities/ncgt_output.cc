#include "../../config.h"

#ifdef USE_NETCDF

#include "ncgt_output.h"
#include <iostream>

using namespace std;

//long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
//		const char* dimension_y, long counter, short update, short rotate_y, double geotop::input::gDoubleNoValue, LONGMATRIX *rc, long **j, long total_pixel){
long ncgt_add_output_var(int ncid, void *m, double time, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
		const char* dimension_y, long counter, short update, short rotate_y, double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j, long total_pixel){


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
	GeoVector<double> tmp_V;

	DOUBLETENSOR *mv;
	DOUBLEMATRIX *mv1;
	DOUBLEMATRIX *m1, *m2;
	GeoMatrix<double> *gm1;  //Alternet to *m1,*m2
	GeoTensor<double> *gm_tens;
	DOUBLEVECTOR *V;
	int r,c,l,i;
	long npoints;
	ncgt_put_double_vs_time(time,dimension_time,counter, ncid,dimension_time);
	//char * function_name="ncgt_add_output_var";
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // TODO (e.g. discharge at the outlet)
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
		if (rotate_y==1){
			ncgt_put_rotate180_y_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time, dimension_x, dimension_y);
		} else {
			ncgt_put_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time, dimension_x, dimension_y);
		}
		break;
	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		gm_tens = new GeoTensor<double>(*((GeoTensor<double> *)m));
		//printf("gm_tens->getDh()=%u,gm_tens->getRh()=%u,gm_tens->getCh()=%u",gm_tens->getDh(),gm_tens->getRh(),gm_tens->getCh());stop_execution();
		tmp_mv.resize(gm_tens->getDh(), gm_tens->getRh(), gm_tens->getCh());
		for (l=1; l<gm_tens->getDh(); l++) {
			for(r=1; r<gm_tens->getRh(); r++) {
				for(c=1; c<gm_tens->getCh(); c++) {
					tmp_mv(l,r,c) = (*gm_tens)[l][r][c];
					}
				}
			}
		tmp_mv.name = gm_tens->name;

		if (rotate_y==1){
			//ncgt_put_rotate180_y_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_rotate180_y_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
		} else {
			//ncgt_put_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
		}
		break;
	case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:// option to print 2D variables in control points
	case NC_GEOTOP_POINT_VAR:
		ncgt_put_doublevector_from_doublematrix_vs_time((DOUBLEMATRIX *)m,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
		break;

	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:// option to print 3D variables in control points
	case NC_GEOTOP_Z_POINT_VAR:
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)m,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		gm_tens = new GeoTensor<double>(*((GeoTensor<double> *)m));
		//printf("gm_tens->getDh()=%u,gm_tens->getRh()=%u,gm_tens->getCh()=%u",gm_tens->getDh(),gm_tens->getRh(),gm_tens->getCh());stop_execution();
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
		ncgt_put_doublematrix_from_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
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
//		mv=new_doubletensor(geotop::common::Variables::Nl,Nr,Nc);
		tmp_mv.resize(geotop::common::Variables::Nl+1, Nr+1, Nc+1);

//		mv->ndl=m2->nrl;
//		for (l=m2->nrl; l<=m2->nrh; l++){
//			for(i=1; i<=total_pixel; i++){
//			//for(i=1; i<=npoints; i++){
//				V->co[i] = m2->co[l][i];
//			}
//			for(r=1; r<= Nr; r++){
//				for (c=1; c<=geotop::common::Variables::Nc; c++){
//		 			if (j[r][c] > 0) {
//						mv->co[l][r][c]=V->co[j[r][c]];
//		 			}else {
//		 				mv->co[l][r][c]=geotop::input::gDoubleNoValue;
//		 			}
//		 		}
//		 	}
//		 }
//		mv->name=m1->name;
		for (l=1; l<gm1->getRows(); l++) {
			for(i=1; i<=total_pixel; i++) {
				tmp_V[i] = (*gm1)[l][i];
				}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=geotop::common::Variables::Nc; c++){
					if (j[r][c] > 0) {
						tmp_mv(l,r,c) = tmp_V[j[r][c]];
					}else {
						tmp_mv(l,r,c) = geotop::input::gDoubleNoValue;
					}
				}
			}
		}
		tmp_mv.name = gm1->name;

		//printf("\nmv->name=%s mv->ndl=%ld",mv->name, mv->ndl);
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		if (rotate_y==1){
			//ncgt_put_rotate180_y_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_rotate180_y_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
		} else {
			//ncgt_put_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
			ncgt_put_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, dimension_x, dimension_y,dimension_z);
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
		//printf("> m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
		printf("> gm1->ncl=%ld, gm1->getCols =%ld,gm1->getRows =%ld,Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",1, gm1->getCols(),gm1->getRows(), Nl, Nr, Nc, total_pixel);
		stop_execution();
		
		tmp_V.resize(npoints+1);
		tmp_mv.resize(geotop::common::Variables::Nl+1, Nr+1, Nc+1);
		
		for (l=1; l<gm1->getRows(); l++) {
			for(i=1; i<=total_pixel; i++) {
				tmp_V[i] = (*gm1)[l][i];
			}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=geotop::common::Variables::Nc; c++){
					if (j[r][c] > 0) {
						tmp_mv(l,r,c) = tmp_V[j[r][c]];
					}else {
						tmp_mv(l,r,c) = geotop::input::gDoubleNoValue;
					}
				}
			}
		}
		tmp_mv.name = gm1->name;
		//printf("name=%s, name1=%s",tmp_mv.name,gm1->name);stop_execution();
		/*
		V = new_doublevector(npoints);
		mv=new_doubletensor(geotop::common::Variables::Nl,Nr,Nc);
		mv->ndl=m2->nrl;
		
		for (l=m2->nrl; l<=m2->nrh; l++){
			for(i=1; i<=total_pixel; i++){
			//for(i=1; i<=npoints; i++){
				V->co[i] = m2->co[l][i];
			}
			for(r=1; r<= Nr; r++){
				for (c=1; c<=geotop::common::Variables::Nc; c++){
					if (j[r][c] > 0) {
						mv->co[l][r][c]=V->co[j[r][c]];
					}else {
						mv->co[l][r][c]=geotop::input::gDoubleNoValue;
					}
				}
			}
		 }
		mv->name=m1->name;
		*/
		//printf("\nmv->name=%s mv->ndl=%ld",mv->name, mv->ndl);
		//ncgt_put_doublematrix_from_doubletensor_vs_time((DOUBLETENSOR *)mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);
		ncgt_put_doublematrix_from_doubletensor_vs_time(tmp_mv,counter, ncid, dimension_time, NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, dimension_z, rc);

		//free_doubletensor(mv);
		//free_doublematrix(m2);
		//free_doublevector(V);
		break;

	case NC_GEOTOP_UNSTRUCT_MAP:
		printf("\ntime=%lf, Hello netcdf maps\n",time);
		m1=(DOUBLEMATRIX*)m;
		m2=new_doublematrix(m1->nrh,m1->nch);
		copy_doublematrix(m1,m2);
		npoints=m1->nch;
		//printf("m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
		mv1=new_doublematrix(Nr,Nc);
		V = new_doublevector(npoints);
		//long npoints=m1->nch;// total_pixel
		for(i=1; i<=total_pixel; i++){
		//for(i=1; i<=npoints; i++){
			V->co[i] = m2->co[1][i];
			}
		for(r=1; r<= Nr; r++){
			for (c=1; c<=geotop::common::Variables::Nc; c++){
				 if (j[r][c] > 0) {
					mv1->co[r][c]=V->co[j[r][c]];
				 }else {
					 mv1->co[r][c]=geotop::input::gDoubleNoValue;
			 	}
			 }
		}
		mv1->name=m1->name;
		printf("\nmv->name=%s",mv1->name);
		//ncgt_put_doublevector_from_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
		if (rotate_y==1){
			ncgt_put_rotate180_y_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time, dimension_x, dimension_y);
		} else {
			ncgt_put_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time, dimension_x, dimension_y);
		}
		free_doublematrix(mv1);
		free_doublematrix(m2);
		free_doublevector(V);
		break;

	case NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT:
		printf("\ntime=%lf, Hello netcdf maps\n",time);
		m1=(DOUBLEMATRIX*)m;
		m2=new_doublematrix(m1->nrh,m1->nch);
		copy_doublematrix(m1,m2);
		npoints=m1->nch;
		//printf("m1->ncl=%ld, m1->nch=%ld, Nl=%ld, Nr=%ld, Nc=%ld, total_pixel=%ld",m1->ncl, m1->nch, Nl, Nr, Nc, total_pixel);stop_execution();
		mv1=new_doublematrix(Nr,Nc);
		V = new_doublevector(npoints);
		//long npoints=m1->nch;// total_pixel
		for(i=1; i<=total_pixel; i++){
		//for(i=1; i<=npoints; i++){
			V->co[i] = m2->co[1][i];
			}
		for(r=1; r<= Nr; r++){
			for (c=1; c<=geotop::common::Variables::Nc; c++){
				 if (j[r][c] > 0) {
					mv1->co[r][c]=V->co[j[r][c]];
				 }else {
					 mv1->co[r][c]=geotop::input::gDoubleNoValue;
				}
			 }
		}
		mv1->name=m1->name;
		printf("\nmv->name=%s",mv1->name);
		ncgt_put_doublevector_from_doublematrix_vs_time((DOUBLEMATRIX *)mv,counter, ncid, dimension_time,  NC_GEOTOP_SUFFIX_FOR_CONTROL_POINT, dimension_x, rc);
		free_doublematrix(mv1);
		free_doublematrix(m2);
		free_doublevector(V);
		break;

	default:
		t_error("incorrect number of dimensions in ncgt_add_output_var");
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
	DOUBLEMATRIX *out2d=NULL;
	DOUBLEMATRIX *out2d0=NULL;
	DOUBLETENSOR *out3d=NULL;
	DOUBLETENSOR *out3d0=NULL;



	//void* m1;// updated matrix
	switch (nlimdim) {
	case NC_GEOTOP_0DIM_VAR: // TODO
		break;
	case NC_GEOTOP_POINT_VAR: // TODO
		break;
	case NC_GEOTOP_2D_MAP:// 2D maps (Y,X)
	case NC_GEOTOP_Z_POINT_VAR:// e.g. point_variable (Z,ID)
	case NC_GEOTOP_2D_MAP_IN_CONTROL_POINT:
		out2d0=(DOUBLEMATRIX*)m0;
		out2d=(DOUBLEMATRIX*)m;
		for (r=out2d0->nrl; r<=out2d0->nrh; r++){
			for (c=out2d0->ncl; c<=out2d0->nch; c++){
				if((out2d0->co[r][c]!=novalue) || ((out2d0->co[r][c]!=out2d0->co[r][c]) && (novalue!=novalue))){
					out2d0->co[r][c]+=out2d->co[r][c]*Dt;
				}
			}
		}
		m0=(void*)out2d0;
		m=(void*)out2d;
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
		out3d0=(DOUBLETENSOR *)m0;
		out3d=(DOUBLETENSOR *)m;
		for (r=out3d0->nrl; r<=out3d0->nrh; r++){
			for (c=out3d0->ncl; c<=out3d0->nch; c++){
				for (l=out3d0->ndl; l<=out3d0->ndh; l++){
					if((out3d0->co[l][r][c]!=novalue) || ((out3d0->co[l][r][c]!=out3d0->co[l][r][c]) && (novalue!=novalue))){
						out3d0->co[l][r][c]+=out3d->co[l][r][c]*Dt;
					}
				}
			}
		}
		m0=(void*)out3d0;
		m=(void*)out3d;
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	default:
		t_error("incorrect number of dimensions in ncgt_var_update");
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

	DOUBLEMATRIX *out2d0=NULL;
	DOUBLETENSOR *out3d0=NULL;
	DOUBLEVECTOR *out1d0=NULL;



	//void* m1;// updated matrix
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
		/*
		out2d0=(DOUBLEMATRIX*)m0;
		for (r=out2d0->nrl; r<=out2d0->nrh; r++){
			for (c=out2d0->ncl; c<=out2d0->nch; c++){
				if((out2d0->co[r][c]!=novalue) || ((out2d0->co[r][c]!=out2d0->co[r][c]) && (novalue!=novalue))){
					out2d0->co[r][c]=0;
				}
			}
		}
		m0=(void*)out2d0;
		*/
		break;

	case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
	case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
		((GeoTensor<double>*)m0)->reset(0, novalue);
		/*
		out3d0=(DOUBLETENSOR*)m0;
		for (r=out3d0->nrl; r<=out3d0->nrh; r++){
			for (c=out3d0->ncl; c<=out3d0->nch; c++){
				for (l=out3d0->ndl; l<=out3d0->ndh; l++){
					if((out3d0->co[l][r][c]!=novalue) || ((out3d0->co[l][r][c]!=out3d0->co[l][r][c]) && (novalue!=novalue))){
						out3d0->co[l][r][c]=0;
					}
				}
			}
		}
		m0=(void*)out3d0;
		*/
		break;
//	printf("\nsono qui1 a=%ld",a);//stop_execution();
	case NC_GEOTOP_UNSTRUCT_MAP:
	case NC_GEOTOP_UNSTRUCT_MAP_IN_CONTROL_POINT:// doublevector
		((GeoVector<double>*)m0)->reset(0, novalue);
		/*
		out1d0=(DOUBLEVECTOR*)m0;
		for (l=out1d0->nl; l<=out1d0->nh; l++){
			if((out1d0->co[l]!=novalue) || ((out1d0->co[l]!=out1d0->co[l]) && (novalue!=novalue))){
				out1d0->co[l]=0;
			}
		}
		m0=(void*)out1d0;
		*/
		break;

	default:
		t_error("incorrect number of dimensions in ncgt_var_set_to_zero");
		break;

	}
	return 0;
}

void* ncgt_new_output_var(const void * m0, const short& nlimdim, const double& novalue, const std::string& suffix, const double& print_flag)
{
   	/* define the temporal counter*/
	/*!
	 * \param m0 - (void *) instantaneous variable (can be doublematrix, doublevector, doubletensor)
	 * \param number_novale - NULL
	 * \param suffix - suffix to be added to the variable_name
	 * \param print_flag - flag on printing option
	 * \description: allocate a new output variable
	 */
	//void* m1;// updated matrix
	void *out=NULL;
	DOUBLEVECTOR *out1d=NULL;
	DOUBLEVECTOR *out1d0=NULL;
	//	DOUBLEMATRIX *out2d=NULL;
	GeoMatrix<double> *out2d=NULL;
	//	DOUBLEMATRIX *out2d0=NULL;
	//	GeoMatrix<double> *out2d0=NULL;
	//DOUBLETENSOR *out3d=NULL;
	//DOUBLETENSOR *out3d0=NULL;
	GeoTensor<double> *out3d = NULL;
	//	GeoTensor<double> *out3d0 = NULL;

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

			//out2d0=(DOUBLEMATRIX *)m0;
			//out2d=new_doublematrix(out2d0->nrh, out2d0->nch);
			//copy_doublematrix(out2d0,out2d);
			//out2d0=out2d;
			//out2d->name=join_strings((char *)out2d0->name,suffix);
			
			return ((void*)out2d);
			break;

		case NC_GEOTOP_3D_MAP:// 3D maps (tensors)
		case NC_GEOTOP_3D_MAP_IN_CONTROL_POINT:
			/*
			out3d0=(DOUBLETENSOR *)m0;
			out3d=new_doubletensor_flexlayer(out3d0->ndl,out3d0->ndh,out3d0->nrh,out3d0->nch);
			out3d->name=join_strings((char *)out3d0->name,suffix);
			copy_doubletensor(out3d0,out3d);
			*/
			out3d = new GeoTensor<double>(*((GeoTensor<double> *)m0));
			out3d->name = out3d->name + suffix;

			//cout << *out3d << endl;
			//cout << "Size of this tensor: " << out3d->getDh() << " x " << out3d->getRh() << " x " << out3d->getCh() << endl;

			//printf("\n out3d->name=%s, check 31/12/2011: DA SISTEMARE\n",out3d->name);
			//(void*)out=out3d;
			return ((void*)out3d);
			break;
	//	printf("\nsono qui1 a=%ld",a);//stop_execution();
		case NC_GEOTOP_UNSTRUCT_MAP:
			out1d0=(DOUBLEVECTOR*)m0;
			out1d=new_doublevector(out1d0->nh);
			copy_doublevector(out1d0,out1d);
			out1d->name=join_strings((char *)out1d0->name,(char*)suffix.c_str());
			break;
		default:
			t_error("incorrect number of dimensions in new_output_var");
			break;

		}
	}
	// initialization of the new variable out
	//ncgt_var_set_to_zero(out,nlimdim, novalue);
	return out;
}

//long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,
//		const char* dimension_y, long counter, short reinitialize, short update, short rotate_y, double geotop::input::gDoubleNoValue, LONGMATRIX *rc, long **j_cont, long total_pixel){
long ncgt_add_output_var_cumtime(int ncid, void *m0, void *m, double time, double computation_time_step, double print_time_step, 
                                 short nlimdim, const char* dimension_time,const char* dimension_z,const char* dimension_x,	const char* dimension_y, 
                                 long counter, short reinitialize, short update, short rotate_y, 
                                 double geotop::input::gDoubleNoValue, GeoMatrix<long>* rc, long **j_cont, long total_pixel){


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
	 * \param rc - (LONGMATRIX) - matrix of the control points
	 * OUTPUT
	 * counter_new: updated counter at which the variable will be written at a successive time
	 */
	long counter_new=counter;
	if(print_time_step>0 && fmod(time,print_time_step)<1.E-5){
		if (m0==NULL) {
			// prints m
			counter_new=ncgt_add_output_var(ncid, m, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter_new,
			update, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
		} else {
			// prints m (instantaneous)
			counter_new=ncgt_add_output_var(ncid, m0, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter,
			NC_GEOTOP_NOUPDATE_COUNTER_TIME, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
			// prints m0 (cumulated) and updates counter
			counter_new=ncgt_add_output_var(ncid, m, time, nlimdim, dimension_time, dimension_z,dimension_x,dimension_y, counter_new,
						update, NC_GEOTOP_ROTATE_Y, NC_GEOTOP_NOVALUE, rc, j_cont, total_pixel);
			// set to zero m0 (cumulated)
			if(reinitialize==1)	ncgt_var_set_to_zero(m0, nlimdim, geotop::input::gDoubleNoValue);
		}
	}else if(print_time_step>0){
		// printing time not reached: updates cumulated variable
		if(reinitialize==1) ncgt_var_update(m, m0, computation_time_step,nlimdim, geotop::input::gDoubleNoValue);
	}
	return counter_new;
}

#endif
