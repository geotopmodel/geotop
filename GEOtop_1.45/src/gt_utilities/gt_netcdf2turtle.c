/*
 * TO BE 
 * 
 */



#ifdef USE_NETCF_ONGOING

#define ERRCODE 2
#define ERROR_MESSAGE(e,n_function,n_ncfunction) {printf("Error in %s() function: %s",n_function,n_ncfunction); printf("\nError: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define GLOBAL_ATTRIBUTE "global_attribute"
#define INIT_VALUE -9998

#include "turtle.h"
#include <netcdf.h>
#include "t_nc_utilities.h"
#include "netcdf2turtle.h"
// README TO GO ON WITH ncgt_NEW_DOUBLEVECTOR



DOUBLEVECTOR *ncgt_new_doublevector(int ncid,const char *varname, const char *dimension){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension (const char *) - name of the dimension
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a doublevector of double[] referred to the variable "varname".
	 */

	const char *function_name="ncgt_new_doublevector";
	size_t lengthp;
	DOUBLEVECTOR *v;
	int status;
	int dimid; /* pointer to the dimenson of the doublevector;*/
	int varid; /* pointer to the variable of the doublevector */



	status=nc_inq_dimid(ncid,dimension,&dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid,&lengthp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	v=new_doublevector((long)lengthp);
	v->name=copy_stringnames(varname);
	status=nc_get_var_double(ncid,varid,&(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_double");

	return v;

}

DOUBLEMATRIX *ncgt_new_doublematrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *
	 *

	*  \param ncid (const char *) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 * \param
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a doublematrix containing the values double[][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_doublematrix";
	size_t lengthp_x,lengthp_y;
	DOUBLEMATRIX *m;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimenson of the doublevector;*/
	int varid; /* pointer to the variable of the doublevector */


	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	m=new_doublematrix((long)lengthp_y,(long)lengthp_x);

	m->name=copy_stringnames(varname);
	status=nc_get_var_double(ncid,varid,&(m->element[m->nrl][m->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_double");



	return m;

}



//EV_S
FLOATVECTOR *ncgt_new_floatvector(int ncid,const char *varname, const char *dimension){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension (const char *) - name of the dimension
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a floatvector of float[] referred to the variable "varname".
	 */

	const char *function_name="ncgt_new_floatvector";
	size_t lengthp;
	FLOATVECTOR *v;
	int status;

	int dimid; /* pointer to the dimenson of the floatvector;*/
	int varid; /* pointer to the variable of the floatvector */


	status=nc_inq_dimid(ncid,dimension,&dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid,&lengthp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	v=new_floatvector((long)lengthp);
	v->name=copy_stringnames(varname);
	status=nc_get_var_float(ncid,varid,&(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_float");



	return v;

}

INTVECTOR *ncgt_new_intvector(int ncid,const char *varname, const char *dimension){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension (const char *) - name of the dimension
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a intvector of int[] referred to the variable "varname".
	 */

	const char *function_name="t_nc_put_intvector";
	size_t lengthp;
	INTVECTOR *v;
	int status;

	int dimid; /* pointer to the dimenson of the intvector;*/
	int varid; /* pointer to the variable of the intvector */




	status=nc_inq_dimid(ncid,dimension,&dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid,&lengthp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	v=new_intvector((long)lengthp);
	v->name=copy_stringnames(varname);
	status=nc_get_var_int(ncid,varid,&(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_int");





	return v;

}

LONGVECTOR *ncgt_new_longvector(int ncid,const char *varname, const char *dimension){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension (const char *) - name of the dimension
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a longvector of long[] referred to the variable "varname".
	 */

	const char *function_name="t_nc_put_longvector";
	size_t lengthp;
	LONGVECTOR *v;
	int status;

	int dimid; /* pointer to the dimenson of the longvector;*/
	int varid; /* pointer to the variable of the longvector */




	status=nc_inq_dimid(ncid,dimension,&dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid,&lengthp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	v=new_longvector((long)lengthp);
	v->name=copy_stringnames(varname);
	status=nc_get_var_long(ncid,varid,&(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_long");





	return v;

}

SHORTVECTOR *ncgt_new_shortvector(int ncid,const char *varname, const char *dimension){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension (const char *) - name of the dimension
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a shortvector of short[] referred to the variable "varname".
	 */

	const char *function_name="t_nc_put_shortvector";
	size_t lengthp;
	SHORTVECTOR *v;
	int status;

	int dimid; /* pointer to the dimenson of the shortvector;*/
	int varid; /* pointer to the variable of the shortvector */




	status=nc_inq_dimid(ncid,dimension,&dimid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid,&lengthp);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	v=new_shortvector((long)lengthp);
	v->name=copy_stringnames(varname);
	status=nc_get_var_short(ncid,varid,&(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_short");





	return v;

}

//CHARVECTOR TBC (unused in Geotop)


FLOATMATRIX *ncgt_new_floatmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a floatmatrix containing the values double[][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_floatmatrix";
	size_t lengthp_x,lengthp_y;
	FLOATMATRIX *m;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimensions of the floatmatrix;*/
	int varid; /* pointer to the variable of the floatmatrix */




	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	m=new_floatmatrix((long)lengthp_y,(long)lengthp_x);

	m->name=copy_stringnames(varname);
	status=nc_get_var_float(ncid,varid,&(m->element[m->nrl][m->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_float");





	return m;

}

SHORTMATRIX *ncgt_new_shortmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a shortmatrix containing the values short[][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_shortmatrix";
	size_t lengthp_x,lengthp_y;
	SHORTMATRIX *m;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimensions of the shortmatrix;*/
	int varid; /* pointer to the variable of the shortmatrix */




	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	m=new_shortmatrix((long)lengthp_y,(long)lengthp_x);

	m->name=copy_stringnames(varname);
	status=nc_get_var_short(ncid,varid,&(m->element[m->nrl][m->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_short");





	return m;

}

INTMATRIX *ncgt_new_intmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a intmatrix containing the values int[][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_intmatrix";
	size_t lengthp_x,lengthp_y;
	INTMATRIX *m;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimensions of the intmatrix;*/
	int varid; /* pointer to the variable of the intmatrix */




	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	m=new_intmatrix((long)lengthp_y,(long)lengthp_x);

	m->name=copy_stringnames(varname);
	status=nc_get_var_int(ncid,varid,&(m->element[m->nrl][m->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_int");





	return m;

}

LONGMATRIX *ncgt_new_longmatrix(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a longmatrix containing the values long[][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_longmatrix";
	size_t lengthp_x,lengthp_y;
	LONGMATRIX *m;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimensions of the longmatrix;*/
	int varid; /* pointer to the variable of the longmatrix */




	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	m=new_longmatrix((long)lengthp_y,(long)lengthp_x);

	m->name=copy_stringnames(varname);
	status=nc_get_var_long(ncid,varid,&(m->element[m->nrl][m->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_long");





	return m;

}

DOUBLETENSOR *ncgt_new_doubletensor(int ncid,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 * \param dimension_z (const char *) - name of the dimension_z
	 *
	 * \author Emanuele Cordano
	 * \date September 2009
	 *
	 *\return a doubletensor containing the values double[][][] referred to the variable "varname".
	 */
//	double *val;
	const char *function_name="ncgt_new_doubletensor";
	size_t lengthp_x,lengthp_y,lengthp_z;
	DOUBLETENSOR *dt;
	int status;

	int dimid_x,dimid_y,dimid_z; /* pointer to the dimensions of the doubletensor;*/
	int varid; /* pointer to the variable of the doubletensor */




	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid x");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid y ");

	status=nc_inq_dimid(ncid,dimension_z,&dimid_z);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid z");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen x");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen y ");

	status=nc_inq_dimlen(ncid,dimid_z,&lengthp_z);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen z ");

	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	dt=new_doubletensor((long)lengthp_x,(long)lengthp_y,(long)lengthp_z);

	dt->name=copy_stringnames(varname);
	status=nc_get_var_double(ncid,varid,&(dt->element[dt->nrl][dt->ncl][dt->ndl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_double");





	return dt;

}
//EV_E



int ncgt_get_doublematrix_level_l(DOUBLEMATRIX *M, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_l) {

	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date  14 October 2008
	 *
	 * \param M (DOUBLEMATRIX *) - doublematrix which is filled (this struct must be allocated and must contain the name of the variable before!)
	 * \param l - (long) - requested layer (0 based)of the variable (M->name)
	 * \param filename (const char *) - name of the NETcdf file
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 * \param dimension_l (const char *) - name of the unlimited dimension (time);
	 *
	 *\brief this functions fills the doublematrix elements with the values of a 3D (time + 2D ) variable M->name contained in the NETCDF file at a certain level (time, index start from 0)
	 *\return 0 in case of success, -1 otherwise
	 *
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 */

	const char *function_name="ncgt_get_doublematrix_level_l";
	size_t lengthp_x,lengthp_y,lengthp_l;
	int status;

	int dimid_x,dimid_y; /* pointer to the dimenson of the doublematrix;*/
	int dimid_l; /* EV pointer to the dimenson of the nc_unlimted variable;*/
	int varid; /* pointer to the variable of the doublematrix */
	long r,c;
	int ndim=3;

	size_t start[ndim],count[ndim];

	if ((M==NULL) || (M->element==NULL)) {
		printf("Error in %s: data struct was not properly allocateted, function returns -1",function_name);
		return -1;
	}



	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	//EV_S
	status=nc_inq_dimid(ncid,dimension_l,&dimid_l);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");
	//EV_E

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	//EV_S
	//status=nc_inq_dimlen(ncid,dimid_y,&lengthp_l); /* length of the time layer dimension written so far */
	status=nc_inq_dimlen(ncid,dimid_l,&lengthp_l); /* length of the time layer dimension written so far */
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");
	//if ((l<0) || (l>(long)lengthp_l)) {
	//The indices are relative to 0,
	if ((l<0) || (l>=(long)lengthp_l)) {
		//EV_E
		printf("Error in %s: variable %s layer %ld is negative or exceeds the maximun leghth %ld, functions returns -1 \n",function_name,M->name,l,(long)lengthp_l);
		return -1;
	}

	status=nc_inq_varid(ncid,M->name,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");
	/* check number of rows and the column */

	//EV_S
	//if (M->nrh!=(long)lengthp_y) {
	if (M->nrh !=(long)lengthp_x) {
		printf("Error in %s: inconstent number of rows %ld (expected) %ld (dimension %s) while reading variable %s, function returns -1 \n",function_name,(long)lengthp_y,M->nrh,dimension_y,M->name);
		return -1;
	}
	//if (M->nch!=(long)lengthp_x) {
	if (M->nch !=(long)lengthp_y) {
		printf("Error in %s: inconstent number of rows %ld (expected) %ld (dimension %s) while reading variable %s, , function returns -1 \n",function_name,(long)lengthp_x,M->nch,dimension_x,M->name);
		return -1;
	}
	//EV_E

	/* initialization of the variable  */
	for(r=M->nrl;r<=M->nrh;r++) {
		for(c=M->ncl;c<=M->nch;c++) {
			M->element[r][c]=INIT_VALUE;
		}
	}
	/* Read the data. Since we know the contents of the file we know
	 * that the data arrays in this program are the correct size to
	 * hold one timestep. */

	count[0]=1;
	count[1]=M->nrh;
	count[2]=M->nch;
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=l;
	start[1]=0;//M->nrl;//EV
	start[2]=0;//M->nrl;//EV

	status=nc_get_vara_double(ncid,varid, start,count, &(M->element[M->nrl][M->ncl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_vara_double");





	/* verify of the correct reading of the variable  */
	for(r=M->nrl;r<=M->nrh;r++) {
		for(c=M->ncl;c<=M->nch;c++) {
			if (M->element[r][c]==INIT_VALUE) {
				printf("Error in %s variable %s element [%ld][%ld] ([%ld][%ld]) was not properly read !!",function_name,M->name,r,c,M->nrh,M->nch);
				return -1;
			}

		}
	}

	return 0;


}

//221009_s
int ncgt_get_doublevector_level_l(DOUBLEVECTOR *v, long l,int ncid,const char *dimension_x, const char *dimension_l){
	/*!
	 *
	 *
	 * \author Enrico Verri
	 * \date  October 2008
	 *
	 * \param v (DOUBLEVECTOR *) - doublevector which is filled (this struct must be allocated and must contain the name of the variable before!)
	 * \param l - (long) - requested layer (0 based)of the variable (v->name)
	 * \param filename (const char *) - name of the NETcdf file
	 * \param dimension_x (const char *) - name of the dimension_x
	  * \param dimension_l (const char *) - name of the unlimited dimension (time);
	 *
	 *\brief this functions fills the doublevector elements with the values of a 2D (time + 1D ) variable v->name contained in the NETCDF file at a certain level (time, index start from 0)
	 *\return 0 in case of success, -1 otherwise
	 *
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 */

	const char *function_name="ncgt_get_doublevector_level_l";
	size_t lengthp_x,lengthp_l;
	int status;

	int dimid_x; /* pointer to the dimenson of the doublevector;*/
	int dimid_l; /* pointer to the dimenson of the nc_unlimted variable;*/
	int varid; /* pointer to the variable of the doublevector */
	long r;
	int ndim=2;

	size_t start[ndim],count[ndim];

	if ((v==NULL) || (v->element==NULL)) {
		printf("Error in %s: vector data struct was not properly allocated, function returns -1",function_name);
		return -1;
	}



	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_l,&dimid_l);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_l,&lengthp_l); /* length of the time layer dimension written so far */
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	//The indices are relative to 0,
	if ((l<0) || (l>=(long)lengthp_l)) {
		//EV_E
		printf("Error in %s: variable %s layer %ld is negative or exceeds the maximun leghth %ld, functions returns -1 \n",function_name,v->name,l,(long)lengthp_l);
		return -1;
	}

	status=nc_inq_varid(ncid,v->name,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");
	/* check number of rows  */
	if (v->nh !=(long)lengthp_x) {
		printf("Error in %s: inconstent number of rows %ld (expected) %ld (dimension %s) while reading variable %s, function returns -1 \n",function_name,(long)lengthp_x,v->nh,dimension_x,v->name);
		return -1;
	}

	/* initialization of the variable  */
	for(r=v->nl;r<=v->nh;r++)  {
			v->element[r]=INIT_VALUE;
	}

	/* Read the data. Since we know the contents of the file we know
	 * that the data arrays in this program are the correct size to
	 * hold one timestep. */
	count[0]=1;
	count[1]=v->nh;
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=l;
	start[1]=0;

	status=nc_get_vara_double(ncid,varid, start,count, &(v->element[v->nl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_vara_double");





	/* verify of the correct reading of the variable  */
	for(r=v->nl;r<=v->nh;r++) {
			if (v->element[r]==INIT_VALUE) {
				printf("Error in %s : variable %s element [%ld] ([%ld]) was not properly read !!",function_name,v->name,r,v->nh);
				return -1;
			}
	}
	return 0;
}

int ncgt_get_doubletensor_level_l(DOUBLETENSOR *t, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l){
	/*!
	 *
	 *
	 * \author Enrico Verri
	 * \date  October 2008
	 *
	 * \param t (DOUBLETENSOR *) - doubletensor which is filled (this struct must be allocated and must contain the name of the variable before!)
	 * \param l - (long) - requested layer (0 based)of the variable (M->name)
	 * \param filename (const char *) - name of the NETcdf file
	 * \param dimension_x (const char *) - name of the dimension_x
	 * \param dimension_y (const char *) - name of the dimension_y
	 * \param dimension_z (const char *) - name of the dimension_z
	 * \param dimension_l (const char *) - name of the unlimited dimension (time);
	 *
	 *\brief this functions fills the doubletensor elements with the values of a 4D (time + 3D ) variable t->name contained in the NETCDF file at a certain level (time, index start from 0)
	 *\return 0 in case of success, -1 otherwise
	 *  //IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  //compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 */

	const char *function_name="ncgt_get_doubletensor_level_l";
	size_t lengthp_x,lengthp_y,lengthp_l,lengthp_z;
	int status;

	int dimid_x,dimid_y,dimid_z; /* pointer to the dimension of the doubletensor;*/
	int dimid_l; /*  pointer to the dimenson of the nc_unlimted variable;*/
	int varid; /* pointer to the variable of the doubletensor */
	long r,c,d;//doubletensor indeces
	int ndim=4;

	size_t start[ndim],count[ndim];

	if ((t==NULL) || (t->element==NULL)) {
		printf("Error in %s: data struct was not properly allocateted, function returns -1",function_name);
		return -1;
	}



	status=nc_inq_dimid(ncid,dimension_x,&dimid_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_y,&dimid_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_z,&dimid_z);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimid(ncid,dimension_l,&dimid_l);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimid");

	status=nc_inq_dimlen(ncid,dimid_x,&lengthp_x);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_y,&lengthp_y);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_z,&lengthp_z);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	status=nc_inq_dimlen(ncid,dimid_l,&lengthp_l); /* length of the time layer dimension written so far */
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_dimlen");

	//The indices are relative to 0,
	if ((l<0) || (l>=(long)lengthp_l)) {
		//EV_E
		printf("Error in %s: variable %s layer %ld is negative or exceeds the maximun leghth %ld, functions returns -1 \n",function_name,t->name,l,(long)lengthp_l);
		return -1;
	}

	status=nc_inq_varid(ncid,t->name,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	/* check number of rows and the column */
	if (t->nrh !=(long)lengthp_x) {
		printf("Error in %s: inconsistent number of rows %ld (expected) %ld (dimension %s) while reading variable %s, function returns -1 \n",function_name,(long)lengthp_y,t->nrh,dimension_y,t->name);
		return -1;
	}
	if (t->nch !=(long)lengthp_y) {
		printf("Error in %s: inconsistent number of columns %ld (expected) %ld (dimension %s) while reading variable %s, , function returns -1 \n",function_name,(long)lengthp_x,t->nch,dimension_x,t->name);
		return -1;
	}
	if (t->ndh !=(long)lengthp_z) {
		printf("Error in %s: inconsistent number of third dimension %ld (expected) %ld (dimension %s) while reading variable %s, , function returns -1 \n",function_name,(long)lengthp_z,t->nch,dimension_z,t->name);
		return -1;
	}

	/* initialization of the variable  */
	for(r=t->nrl;r<=t->nrh;r++) {
		for(c=t->ncl;c<=t->nch;c++) {
			for(d=t->ndl;d<=t->ndh;d++) {
				t->element[r][c][d]=INIT_VALUE;
			}
		}
	}

	/* Read the data. Since we know the contents of the file we know
	 * that the data arrays in this program are the correct size to
	 * hold one timestep. */
	count[0]=1;
	count[1]=t->nrh;
	count[2]=t->nch;
	count[3]=t->ndh;
	//The indices are relative to 0, so for example,the first data value of a variable would have index (0, 0, ... , 0).
	start[0]=l;
	start[1]=0;
	start[2]=0;
	start[3]=0;

	status=nc_get_vara_double(ncid,varid, start,count, &(t->element[t->nrl][t->ncl][t->ndl]));
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_vara_double");





	/* verify of the correct reading of the variable  */
	for(r=t->nrl;r<=t->nrh;r++){
		for(c=t->ncl;c<=t->nch;c++){
			for(d=t->ndl;d<=t->ndh;d++){
				if (t->element[r][c][d]==INIT_VALUE){
					printf("Error in %s variable %s element [%ld][%ld][%ld] ([%ld][%ld][%ld]) was not properly read !!",function_name,t->name,r,c,d,t->nrh,t->nch,t->ndh);
					return -1;
				}
			}
		}
	}
	return 0;
}
//221009_e

//201009_s
//SCALAR TYPES
//
//NC_BYTE
signed char ncgt_get_byte(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a byte referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="ncgt_get_byte";//
	signed char bvar;//
	int status;

	int varid; /* pointer to the variable of the byte */



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_schar(ncid, varid, &bvar);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_schar");//




	return bvar;

}

// NC_INT
int nc_get_int(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a int32 referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="nc_get_int";//
	int ivar;//
	int status;

	int varid; /* pointer to the variable of the int *///



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_int(ncid, varid, &ivar);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_int");//




	return ivar;//

}

// NC_SHORT
short nc_get_short(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a short16 referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="nc_get_short";//
	short svar;//
	int status;

	int varid; /* pointer to the variable of the short *///



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_short(ncid, varid, &svar);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_short");//




	return svar;//

}

// NC_FLOAT
float nc_get_float(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a float32 referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="nc_get_float";//
	float fvar;//
	int status;

	int varid=0; /* pointer to the variable of the float *///



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_float(ncid, varid, &fvar);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_float");//




	return fvar;//

}
// NC_DOUBLE
double nc_get_double(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a double64 referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="nc_get_double";//
	double d64var = 0.0;//
	int status;

	int varid=0; /* pointer to the variable of the double *///



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_double(ncid, varid, &d64var);//
	//status = nc_get_var1_double(ncid, varid,0, &d64var);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_double");//




	return d64var;//

}

// NC_LONG
long nc_get_long(int ncid,const char *varname)
{
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param varname  (const char *) - name of the variable
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\return a long64 referred to the value of the netcdf variable "varname".//
	 */

	const char *function_name="nc_get_long";//
	long lvar;//
	int status;

	int varid; /* pointer to the variable of the long *///



	status=nc_inq_varid(ncid,varname,&varid);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_inq_varid");

	status = nc_get_var_long(ncid, varid, &lvar);//
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_var_long");//




	return lvar;//

}
//201009_e
//131109_s
DOUBLETENSOR *ncgt_new_doubletensor_rotate180y(int ncid,const char *varname, const char *dimension_x,const char *dimension_y,const char *dimension_z){
	/*!
	 *\param fileneme - (char *) name of the NetCDF file
	 *\param varname  name of the variable
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *\param dimension_z - (char *) name of the z dimension
	 *
	 *\brief This function read and rotates the variable contained in a netcdf
	 *\brief using the functions rotate180_y_doubletensor() and ncgt_new_doubletensor() .
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 * \author Enrico Verri,Emanuele Cordano
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.

	 *
	 */

	DOUBLETENSOR *dt;

	char *function_name="ncgt_new_doubletensor_rotate180y";
	int res = 0;

	dt=ncgt_new_doubletensor(filename,varname,dimension_x,dimension_y,dimension_z);
	res=rotate180_y_doubletensor(dt);
	if (res!=0) printf("Error in %s: rotate180_y_doubletensor() did not work correctly  \n",function_name);

	return dt;
}
DOUBLEMATRIX *ncgt_new_doublematrix_rotate180y(int ncid,const char *varname, const char *dimension_x,const char *dimension_y){
	/*!
	 *\param fileneme - (char *) name of the NetCDF file
	 *\param varname  name of the variable
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *
	 *\brief This function read and rotates the variable contained in a netcdf
	 *\brief using the functions rotate180_y_doubletensor() and ncgt_get_doubletensor() .
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 * \author Enrico Verri,Emanuele Cordano
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.

	 *
	 */

	DOUBLEMATRIX *dm;

	char *function_name="ncgt_new_doubletensor_rotate180y";
	int res = 0;

	dm=ncgt_new_doublematrix(filename,varname,dimension_x,dimension_y);
	res=rotate180_y_doublematrix(dm);
	if (res!=0) printf("Error in %s: rotate180_y_doublematrix() did not work correctly  \n",function_name);

	return dm;

}
int ncgt_get_doubletensor_rotate180y_level_l(DOUBLETENSOR *t, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_z, const char *dimension_l){
	/*!
	 *\param t - (DOUBLETENSOR *) variable from the NetCDF
	 *\param fileneme - (char *) name of the NetCDF file
	 *\param l        - (long) number of the level (0 based) at which the xyz map is printed
	 *\param dimension_t - (char *) name of the t dimension (number of times: UNLIMITED)
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *\param dimension_z - (char *) name of the z dimension
	 *
	 *\brief This function read and rotates the variable contained in a netcdf
	 *\brief referred to a particular time step 3D + time variable
	 *\brief using the functions rotate180_y_doubletensor() and ncgt_get_doubletensor_level_l() .
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 * \author Enrico Verri,Emanuele Cordano
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.

	 *
	 */

	int res=0;

	char *function_name="ncgt_get_doubletensor_rotate180y_level_l";

	res=ncgt_get_doubletensor_level_l(t,l,filename,dimension_x,dimension_y,dimension_z,dimension_l);
	if (res!=0) printf("Error in %s: ncgt_get_doubletensor_level_l() did not work correctly  \n",function_name);
	res=rotate180_y_doubletensor(t);
	if (res!=0) printf("Error in %s: rotate180_y_doubletensor() did not work correctly  \n",function_name);

	return res;


}
int ncgt_get_doublematrix_rotate180y_level_l(DOUBLEMATRIX *M, long l,int ncid,const char *dimension_x, const char *dimension_y, const char *dimension_l){
	/*!
	 *\param M - (DOUBLEMATRIX *) variable from the NetCDF
	 *\param fileneme - (char *) name of the NetCDF file
	 *\param l        - (long) number of the level (0 based) at which the xy map is printed
	 *\param dimension_l - (char *) name of the t dimension (number of times: UNLIMITED)
	 *\param dimension_x - (char *) name of the x dimension (number of column)
	 *\param dimension_y - (char *) name of the y dimension (number of row)
	 *
	 *\brief This function read and rotates the variable contained in a netcdf
	 *\brief referred to a particular time step 2D + time variable
	 *\brief using the functions rotate180_y_doublematrix() and ncgt_get_doublematrix_level_l() .
	 *  IMPORTANT: only one NC_UNLIMITED dimension allowed
	 *  compile netcdf_lib with --enable-netcdf-4 option tu Turn on netCDF-4 features.
	 *
	 * \author Enrico Verri,Emanuele Cordano
	 * \date October 2009
	 *
	 * \return 0 if exit is ok, otherwise an error message.

	 *
	 */

	int res=0;

	char *function_name="ncgt_get_doublematrix_rotate180y_level_l";

	res=ncgt_get_doublematrix_level_l(M,l,filename,dimension_x,dimension_y,dimension_l);
	if (res!=0) printf("Error in %s: ncgt_get_doublematrix_level_l() did not work correctly  \n",function_name);
	res=rotate180_y_doublematrix(M);
	if (res!=0) printf("Error in %s: rotate180_y_doublematrix() did not work correctly  \n",function_name);

	return res;

}
void nc_get_global_attr_lat_lon_min_max(int ncid,double *long_min,double *long_max,double *lat_min,double *lat_max){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\load min max lat & long netcdf global attributes in local variables.//
	 */

	const char *function_name="nc_get_global_attr_lat_lon_min_max";//
	int status;

	//int varid; /* pointer to the variable of the byte */



	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_LONG_MIN, long_min);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//
	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_LONG_MAX, long_max);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//
	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_LAT_MIN, lat_min);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//
	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_LAT_MAX, lat_max);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//





}

//131109_e
void nc_get_global_attr_missing_value(int ncid,double *missing_value){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param missing_value (double*)
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\load missing_value netcdf global attribute in local variable.//
	 *\To handle NaN special value see dumplib.c and isnan.h in ncdump netcdf sample
	 */

	const char *function_name="nc_get_global_attr_missing_value";//
	int status;




	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_MISSING_VALUE, missing_value);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//




}

void nc_get_global_attr_resolution(int ncid,double *resolution){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param resolution (double*)
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\load map resolution netcdf global attribute in local variable.//
	 */

	const char *function_name="nc_get_global_attr_resolution";//
	int status;




	status = nc_get_att_double(ncid, NC_GLOBAL, GLOB_ATTR_MAP_RESOLUTION, resolution);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//





}

void nc_get_global_value(int ncid,char *attr_name,double *attr_value){
	/*!
	 *
	 *

	*  \param ncid (int) - pointer to the netcdf file
	 * \param resolution (double*)
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 *
	 *\load map resolution netcdf global attribute in local variable.//
	 */

	const char *function_name="nc_get_global_value";//
	int status;




	status = nc_get_att_double(ncid, NC_GLOBAL, attr_name, attr_value);
	if (status!=NC_NOERR) ERROR_MESSAGE(status,function_name,"nc_get_att_double");//





}

#endif
