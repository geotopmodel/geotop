
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "struct.geotop.h"
#include "constants.h"
#include "meteodata.h"
#include "times.h"
#include "meteo.h"
#include "../libraries/ascii/rw_maps.h"
#include "../libraries/ascii/tabs.h"

extern long number_absent, number_novalue;
extern char *string_novalue;


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_interp(short flag, short update, long *istart, double *out, double **data, long nlines, long ncols, long col_date, double tbeg, double tend)

{
	
	short abeg=0, aend=0;
	long ibeg, iend, i, j;
	double t;
	
	ibeg = find_line_data(flag, tbeg, *istart, data, col_date, nlines, &abeg);
	iend = find_line_data(flag, tend, ibeg, data, col_date, nlines, &aend);
	
	//printf("istart:%ld ibeg:%ld iend:%ld %f %f %ld %ld\n",*istart,ibeg,iend,tbeg,tend,abeg,aend);
	
	for (i=0; i<ncols; i++) {
		
		if ( (long)data[0][i] == number_absent) {
			
			out[i] = (double)number_absent;
			
			//printf("->i:%ld out:%f\n\n",i,out[i]);
			
		}else if (abeg == 1 && aend == 1) {
			
			out[i] = 0.;
			
			if ( (long)out[i] != number_novalue) {
				out[i] += integrate_meas(flag, tbeg, ibeg+1, data, i, col_date);
				//printf("0:tbeg:%f %f\n",tbeg,integrate_meas(flag, tbeg, ibeg+1, data, i, col_date));
			}
			
			if ( (long)out[i] != number_novalue) {
				out[i] -= integrate_meas(flag, tend, iend+1, data, i, col_date);
				//printf("end:tend:%f %f\n",tend,-integrate_meas(flag, tend, iend+1, data, i, col_date));
			}
			
			for (j=ibeg+1; j<=iend; j++) {
				if ( (long)out[i] != number_novalue) {
					t = time_in_JDfrom0(flag, j, col_date, data);
					out[i] += integrate_meas(flag, t, j+1, data, i, col_date);
					//printf("j:%ld:t:%f %f\n",j,t,integrate_meas(flag, t, j+1, data, i, col_date));
				}
			}
			
			if ( (long)out[i] != number_novalue) {
				out[i] /= (tend	- tbeg);
			}
			
			//printf("i:%ld out:%f\n\n",i,out[i]);
			
		}else {
			
			out[i] = (double)number_novalue;
			
			//printf("<-i:%ld out:%f\n\n",i,out[i]);
			
		}
	}
	
	if (update>0) *istart = ibeg;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void meteo_interp2(short flag, long *istart, double *out, double **data, long nlines, long ncols, long col_date, double tbeg, double tend)

{
	
	short abeg;
	long ibeg, i;
	
	ibeg = find_line_data(flag, tbeg, *istart, data, col_date, nlines, &abeg);
	if (abeg == 3){
		ibeg = nlines;
		abeg = 1;
	}
	
	for (i=0; i<ncols; i++) {
		
		if ( (long)data[0][i] == number_absent) {
			
			out[i] = (double)number_absent;
			
		}else if (abeg == 1 && (long)data[ibeg][i] != number_novalue) {
			
			out[i] = data[ibeg][i];
			
		}else {
			
			out[i] = (double)number_novalue;
			
		}
	}
	
	
	*istart = ibeg;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fixing_dates(long imeteo, double **data, double ST, double STstat, long nlines, long date12col, long JDfrom0col){
	
	long i;
	double dateold=0.0;
	
	if ( (long)data[0][JDfrom0col] == number_absent) {
		
		for (i=0; i<nlines; i++) {
			//converting in JDfrom0
			data[i][JDfrom0col] = convert_dateeur12_JDfrom0(data[i][date12col]);
			//setting ST 
			data[i][JDfrom0col] += (ST - STstat) / 24.;
			
			if(data[i][JDfrom0col]<=dateold){
				printf("Error: at line %ld meteo file %ld time is equal or earlier than previous line %f %f\n",i+1,imeteo,dateold,data[i][JDfrom0col]);
				stop_execution();
				t_error("Not Possible To Continue");
			} 		
			dateold = data[i][JDfrom0col];
		}
		return 1;
	}else {
		return 0;
	}	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_wind_speed(double **data, long nlines, long Wspeed, long Wdir, long Wx, long Wy, char *HeaderWx, char *HeaderWy){
	
	long i;
	
	//if the columns Wspeed and Wdir are present, and the columns Wx and Wy are not present
	if ( (long)data[0][Wspeed] != number_absent && (long)data[0][Wdir] != number_absent && ( (long)data[0][Wx] == number_absent || (long)data[0][Wy] == number_absent ) ) {
		for (i=0; i<nlines; i++) {
			if ( (long)data[i][Wspeed] != number_novalue && (long)data[i][Wdir] != number_novalue ) {
				data[i][Wx] = -data[i][Wspeed] * sin(data[i][Wdir] * Pi / 180.);
				data[i][Wy] = -data[i][Wspeed] * cos(data[i][Wdir] * Pi / 180.);
			}else {
				data[i][Wx] = (double)number_novalue;
				data[i][Wy] = (double)number_novalue;
			}
			
			if(strcmp(HeaderWx,string_novalue)!=0 && strcmp(HeaderWy,string_novalue)!=0){
				data[i][Wspeed] = (double)number_absent;
				data[i][Wdir] = (double)number_absent;	
			}			
		}
		if(strcmp(HeaderWx,string_novalue)!=0 && strcmp(HeaderWy,string_novalue)!=0){
			return 1;
		}else {
			return 0;
		}
	}else {
		return 0;
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short fill_Tdew(long imeteo, DOUBLEVECTOR *Z, double **data, long nlines, long RH, long Tair, long Tairdew, char *HeaderTdew, double RHmin){
	
	long i;
	
	if ( (long)data[0][RH] != number_absent && (long)data[0][Tair] != number_absent && (long)data[0][Tairdew] == number_absent ) {
		for (i=0; i<nlines; i++) {
			if ( (long)data[i][RH] != number_novalue && (long)data[i][Tair] != number_novalue ) {
				data[i][Tairdew] = Tdew(data[i][Tair], Fmax(RHmin, data[i][RH])/100., Z->co[imeteo]);
			}else {
				data[i][Tairdew] = (double)number_novalue;
			}
			
			if (strcmp(HeaderTdew, string_novalue) != 0) {
				data[i][RH] = (double)number_absent;
			}
		}
		if (strcmp(HeaderTdew, string_novalue) != 0) {
			return 1;
		}else {
			return 0;
		}
	}else {
		return 0;
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rewrite_meteo_files(double **meteo, long meteolines, char **header, char *name, short added_JD, short added_wind, short added_cloudiness, short added_Tdew){
	
	char *newname;
	short first_column, write;
	long i, n, d, m, y, h, mi;
	FILE *f;
	
	newname = join_strings(name, ".old");
	f=fopen(newname,"r");
	if (f == NULL){
		write = 1;
	}else {
		fclose(f);
		write = 0;
	}
	
	if (added_cloudiness == 0 && added_wind == 0 && added_JD == 0 && added_Tdew == 0) write = 0;
	
	if (write == 1) {
		rename(name, newname);
		
		f = fopen(name, "w");
		
		first_column = 1;
		for (i=0; i<nmet; i++) {
			if ( (long)meteo[0][i] != number_absent && strcmp(header[i], string_novalue) != 0){
				if (first_column == 0) {
					fprintf(f,",");
				}else {
					first_column = 0;
				}					
				fprintf(f, "%s", header[i]);
			}
		}
		fprintf(f, "\n");
		
		for (n=0; n<meteolines; n++) {
			first_column = 1;
			for (i=0; i<nmet; i++) {
				if ( (long)meteo[0][i] != number_absent && strcmp(header[i], string_novalue) != 0){
					if (first_column == 0) {
						fprintf(f,",");
					}else {
						first_column = 0;
					}
					if(i == iDate12) {
						convert_dateeur12_daymonthyearhourmin(meteo[n][i], &d, &m, &y, &h, &mi);
						fprintf(f, "%02.0f/%02.0f/%04.0f %02.0f:%02.0f",(float)d,(float)m,(float)y,(float)h,(float)mi);
					}else{
						fprintf(f, "%f", meteo[n][i]);
					}
				}
			}
			fprintf(f, "\n");
		}
		
		fclose(f);
	}
	
	free(newname);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double integrate_meas(short flag, double t, long i, double **data, long col, long col_date){
	
	double t0, t1, value, res;
	
	if( (long)data[i  ][col] != number_novalue && (long)data[i  ][col] != number_absent && 
	   (long)data[i-1][col] != number_novalue && (long)data[i-1][col] != number_absent ) {
		
		t0 = time_in_JDfrom0(flag, i-1, col_date, data);
		t1 = time_in_JDfrom0(flag, i  , col_date, data);
		
		if(fabs(t0-t1) < 1.E-5){
			printf("There are 2 consecutive line in a meteo file with same date: t0:%f(%ld) t1:%f(%ld) equal",t0,i-1,t1,i);
			stop_execution();
			t_error("Not Possible To Continue");
		}
		
		value = ( (t - t0) * data[i][col] + (t1 - t) * data[i-1][col] )/(t1 - t0);
		//printf("Integrate: t:%f t0:%f t1:%f v0:%f v1:%f v:%f\n",t,t0,t1,data[i-1][col],data[i][col],value);
		
		//area of trapezium 
		res = 0.5 * ( value + data[i][col] ) * (t1 - t);	
		
	}else {
		
		res = (double)number_novalue;
		
	}
	
	return(res);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long find_line_data(short flag, double t, long ibeg, double **data, long col_date, long nlines, short *a){
	
	long i;
	double t0, t1;
	
	i = ibeg;
	
	do{
		*a = 0;
		
		if (i+1 <= nlines-1) {
			
			t0 = time_in_JDfrom0(flag, i  , col_date, data);
			t1 = time_in_JDfrom0(flag, i+1, col_date, data);
			
			if (t0 <= t && t1 >= t){
				*a = 1;
			}else if (t0 > t){
				*a = 2;
			}
			
		}
		
		if (i+1 >= nlines-1 && *a == 0) *a = 3;
		
		if (*a == 0) i++;
		
	}while (*a == 0);
	
	return(i);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

// All the times in the following subroutines are in Julian days from year 0

// flag 0: dates read in DDMMYYYYhhmm and converted in JDfrom0
// flag 1: dates in JDfrom0

double time_in_JDfrom0(short flag, long i, long col, double **data){
	
	double t;
	
	if (flag == 0) {
		t = convert_dateeur12_JDfrom0(data[i][col]);
	}else {
		t = data[i][col];
	}
	
	return(t);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

long find_station(long metvar, long nstat, double **var){
	
	long i=0;
	
	while ( (long)var[i][metvar] == number_absent && i < nstat-1 ) {
		i++;
	}
	
	return i;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double **read_horizon(short a, long i, char *name, char **ColDescr, long *num_lines, FILE *flog){
	
	FILE *f;
	long j;
	double **hor;
	char *temp;
	short fileyes;
	
	//check is the file exists
	if(strcmp(name , string_novalue) != 0){
		fileyes=1;
	}else{
		fileyes=-1;
	}
	if(fileyes==1){
		temp=namefile_i(name,i);
		f=fopen(temp,"r");
		if(f==NULL){
			fileyes=0;
		}else{
			fclose(f);
		}
		free(temp);
	}
		
	//different cases
	if (fileyes == -1) {
		
		if(a==0){
			fprintf(flog,"\nWarning:: No horizon file found for point type #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			printf("\nWarning:: No horizon file found for point type #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		}else if (a==1){
			fprintf(flog,"\nWarning:: No horizon file found for meteo station #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			printf("\nWarning:: No horizon file found for meteo station #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		}
		
		*num_lines = 4;
		hor = (double**)malloc((*num_lines)*sizeof(double*));
		for( j=0; j<(*num_lines); j++){
			hor[j] = (double*)malloc(2*sizeof(double));
			hor[j][0] = 45.0+j*90.0;
			hor[j][1] = 0.0;
		}
		
	}else if (fileyes == 0) {
		
		if(a==0){
			fprintf(flog,"\nWarning:: No horizon file found for point type #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			printf("\nWarning:: No horizon file found for point type #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		}else if (a==1){
			fprintf(flog,"\nWarning:: No horizon file found for meteo station #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
			printf("\nWarning:: No horizon file found for meteo station #%ld. In this case the horizon will be considered always not obscured, i.e. shadow=FALSE\n",i);
		}
		
		temp=namefile_i(name,i);
		f=fopen(temp,"w");
		fprintf(f,"! Horizon file for met station or point #%ld \n",i);
		fprintf(f,"! All measures in degrees\n");
		fprintf(f,"\n");
		fprintf(f,"%s,%s\n",ColDescr[0],ColDescr[1]);
		
		*num_lines = 4;
		hor = (double**)malloc((*num_lines)*sizeof(double*));
		for( j=0; j<(*num_lines); j++){
			hor[j] = (double*)malloc(2*sizeof(double));
			hor[j][0] = 45.0+j*90.0;
			hor[j][1] = 0.0;
			fprintf(f,"%f,%f\n",hor[j][0],hor[j][1]);
		}
		
		fclose(f);
		free(temp);
		
	} else if(fileyes==1){
		
		if(a==0){
			fprintf(flog,"\nHorizon file FOUND for point type #%ld\n",i);
			printf("\nHorizon file FOUND for point type #%ld\n",i);
		}else if (a==1) {
			fprintf(flog,"\nHorizon file FOUND for meteo station #%ld\n",i);
			printf("\nHorizon file FOUND for meteo station #%ld\n",i);
		}
		
		temp = namefile_i(name,i);
		hor = read_txt_matrix(temp, 33, 44, ColDescr, 2, num_lines, flog);
		free(temp);
		
		if ( (long)hor[0][0] == number_absent || (long)hor[0][1] == number_absent) {
			printf("Error:: In the file %s the columns %s and/or %s are missing\n",temp,ColDescr[0],ColDescr[1]);
			stop_execution();
			t_error("Not Possible To Continue (3) ");
		}
		
	}
	
	return(hor);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
