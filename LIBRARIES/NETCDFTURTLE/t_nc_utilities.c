
/* Turtle_NetCdf CONTAINS FUNCTIONS TO INPORT/EXPORT FUIDTURLE STRUCT OF DATA IN NETcdf FILES AS VARIABLES
Turtle_NetCdf Version 0.9375 KMackenzie

file t_nc_utilities.c

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of numerioc_solver.
	 Turtle_NetCdf is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

     Turtle_NetCdf is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 *
 * \file t_nc_utilities.c
 * \data 15 September 2009
 * \author Emanuele Cordano
 */

#include "turtle.h"
#include "t_nc_utilities.h"
#include "turtle2netcdf.h" //061009

char *copy_stringnames(const char *origin){
/*!
 * \author Emanuele Cordano
 * \date April 2009
 *
 *
 */
char *dest;

//dest=origin;
dest=(char *)malloc((size_t)(strlen(origin)+4)*sizeof(char));
strcpy(dest,origin);
return dest;
}

//061009_s
//http://cf-pcmdi.llnl.gov/documents/cf-standard-names/standard-name-table/12/cf-standard-name-table.html
void nc_add_cf_convention_header(const char *filename,const char *title, const char *institution,const char *source,const char *history,const char *references,const char *comment,const char *conventions){
	/*!
	 *
	 * netcdf CF standard provides some basic discovery metadata in global attributes
	 * \param filename (const char *) - name of the nc file
	 * \param title (const char *) - Whats in the file (A succinct description of what is in the dataset)
	 * \param institution (const char *) - Where it was produced
	 * \param source (const char *) - How it was produced e.g. model version, instrument type
	 * \param history (const char *) -Provides an audit trail for modifications to the original data. Well-behaved generic netCDF filters will
	 *								  automatically append their name and the parameters with which they were invoked to the global history attribute
	 *								  of an input netCDF file. We recommend that each line begin with a timestamp indicating the date and time of
	 *								  day that the program was executed
	 * \param references (const char *) - Published or web-based references that describe the data or methods used to produce it.
	 * \param comment (const char *) - Miscellaneous information about the data or methods used to produce it.
	 * \param conventions (const char *) - Name of the conventions followed by the dataset (I.E. 1.0 .. 1.4)
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 */

	//DO NOT USE ":" at the end of attribute_name parameters
	t_nc_put_var_textattributes(filename,"global_attribute", "Conventions",conventions);
	t_nc_put_var_textattributes(filename,"global_attribute", "title",title);
	t_nc_put_var_textattributes(filename,"global_attribute", "institution",institution);
	t_nc_put_var_textattributes(filename,"global_attribute", "source",source);
	t_nc_put_var_textattributes(filename,"global_attribute", "history",history);
	t_nc_put_var_textattributes(filename,"global_attribute", "references",references);
	t_nc_put_var_textattributes(filename,"global_attribute", "comment",comment);


}
//
void nc_add_cf_convention_var_attributes(const char *filename,const char *varname, const char *units,const char *standard_name,const char *long_name,const char *ancillary_variables,
										 const char *fillvalue,const char *missing_value){//,const char *valid_max,const char *valid_min,const char *valid_range){
	/*!
	 *
	 * netcdf CF standard provides some basic local attributes (mandatory)
	 * \param filename (const char *) - name of the nc file
	 * \param varname (const char *) - name of the nc variable
	 * \param units (const char *) - is mandatory for all variables containing data other than dimensionless
	 * 		numbers. The units do not identify the physical quantity. units can be
	 *  		udunits strings e.g. 1, degC, Pa, mbar, W m-2, kg/m2/s, mm day^-1
	 *		or COARDS specials layer, level, sigma level. Udunits doesnt support ppm, psu, dB, Sv.
	 * \param standard_name (const char *) - identifies the quantity. Units must be consistent with standard
	 * 		name and any statistical processing e.g. variance.
	 *		standard_name does not include coordinate or processing information.
	 * \param long_name (const char *) -is not standardised.
	 *		ancillary _ariables (const char *) -is a pointer to variables providing metadata about
	 * 		the individual data values e.g. standard error or data quality information.
	 * Numeric data variables may have FillValue, missing value, valid max,
	 * valid min, valid range. CF deprecates missing value.
	 * Variables containing flag values need flag values and flag meanings to
	 * make them self-describing.
	 *
	 * \author Enrico Verri
	 * \date October 2009
	 */

	t_nc_put_var_textattributes(filename,varname, "units",units);
	t_nc_put_var_textattributes(filename,varname, "standard_name",standard_name);
	t_nc_put_var_textattributes(filename,varname, "long_name",long_name);
	t_nc_put_var_textattributes(filename,varname, "ancillary_variables",ancillary_variables);
	t_nc_put_var_textattributes(filename,varname, "fillvalue",fillvalue);
	/*t_nc_put_var_textattributes(filename,varname, "valid_max",valid_max);
	t_nc_put_var_textattributes(filename,varname, "valid_min",valid_min);
	t_nc_put_var_textattributes(filename,varname, "valid_range",valid_range);*/

}

//061009_e


int rotate180_y_doublematrix(DOUBLEMATRIX *M) {
/*
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLEMATRIX *M) matrix to rotate
 */
	long r,c;
	double val=0.0;


	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_doubletensor(DOUBLETENSOR *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLETENSOR *M) matrix to rotate
 */
	long r,c,l;
	double val=0.0;

	long r_max=(M->nrh-M->nrl+1)/2;
//	printf("r_max=%ld %ld %ld \n",r_max,M->nrh,M->nrl);

	for(l=M->ndl;l<=M->ndh;l++) {
		for(r=M->nrl;r<=r_max;r++) {
			for (c=M->ncl;c<=M->nch;c++) {
				val=M->element[l][r][c];
				M->element[l][r][c]=M->element[l][M->nrh-r+1][c];
				M->element[l][M->nrh-r+1][c]=val;
			}
		}
	}
	return 0;
}

int rotate180_y_floatmatrix(FLOATMATRIX *M) {
/*
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (DOUBLEMATRIX *M) matrix to rotate
 */
	long r,c;
	float val=0.0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}



int rotate180_y_longmatrix(LONGMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (LONGMATRIX *M) matrix to rotate
 */
	long r,c;
	long val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_intmatrix(INTMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (INTMATRIX *M) matrix to rotate
 */
	long r,c;
	int val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int rotate180_y_shortmatrix(SHORTMATRIX *M) {
/*!
 * \author Emanuele Cordano
 * \date October 2009
 *
 * \param (SHORTMATRIX *M) matrix to rotate
 */
	long r,c;
	short val=0;

	long r_max=(M->nrh-M->nrl+1)/2;

	for(r=M->nrl;r<=r_max;r++) {
		for (c=M->ncl;c<=M->nch;c++) {
			val=M->element[r][c];
			M->element[r][c]=M->element[M->nrh-r+1][c];
			M->element[M->nrh-r+1][c]=val;

		}
	}

	return 0;
}

int invert_order_doublevector(DOUBLEVECTOR *v) {
/*!
 *
 * \author Emanuele Cordano
 * \date October 2009
 *
 *\param (DOUBLEVECTOR *) the vector which is inverted
 */

	long j;
	long j_max=(v->nh-v->nl+1)/2;
	double val=0.0;
	for (j=v->nl;j<=j_max;j++) {
		val=v->element[j];
		v->element[j]=v->element[v->nh-j+1];
		v->element[v->nh-j+1]=val;
	}

	return 0;
}

int invert_order_longvector(LONGVECTOR *v) {
/*!
 *
 * \author Emanuele Cordano
 * \date October 2009
 *
 *\param (LONGVECTOR *) the vector which is inverted
 */

	long j;
	long j_max=(v->nh-v->nl+1)/2;
	long val=0;
	for (j=v->nl;j<=j_max;j++) {
		val=v->element[j];
		v->element[j]=v->element[v->nh-j+1];
		v->element[v->nh-j+1]=val;
	}

	return 0;
}

