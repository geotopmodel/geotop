
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
#include "clouds.h"
#include "constants.h"
#include "meteo.h"
#include "radiation.h"
#include "times.h"

#define ndivday 3
#define filecloud "clouds.txt"

extern long number_novalue, number_absent;
extern char *WORKING_DIRECTORY;

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//lat and lon in [deg]
short fill_meteo_data_with_cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat, 
									 double lon, double ST, double Z, double sky, double A){
	
	double *cloudtrans;
	long n;
	
	//if there are radiation data, and no cloudiness
	if ( ( (long)meteo[0][iSW] != number_absent || ( (long)meteo[0][iSWb] != number_absent && (long)meteo[0][iSWd] != number_absent) ) && 
	     (long)meteo[0][itauC] == number_absent && (long)meteo[0][iC] == number_absent ){
			
		cloudtrans = (double*)malloc(meteolines*sizeof(double));
		cloudiness(meteo, meteolines, horizon, horizonlines, lat*Pi/180., lon*Pi/180., ST, Z, sky, A, cloudtrans);
				
		for (n=0; n<meteolines; n++) {
			meteo[n][itauC] = cloudtrans[n];
		}
			
		free(cloudtrans);
		return 1;	
	}else {
		return 0;
	}


}
	
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

void cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat, double lon, 
				double ST, double Z, double sky, double A, double *cloudtrans){
	
	double tc, tc0;
	double E0, Et, Delta, alpha, direction;
	long n00, n0, n1, k, n;
	long *ndiv;
	FILE *f;
	char *temp;
	
	//file header
	temp = join_strings(WORKING_DIRECTORY,filecloud);
	f = fopen(temp,"w");
	fprintf(f,"Date,SolarHeight[deg],SolarAzimuth[deg],SinSolarHeight,SWinMeasured[W/m2],SWinClearSky[W/m2],AtmTransmissivity,CloudTransmissivity\n");
	fclose(f);	
	free(temp);
	
	//average cloudiness in ndivday daily intervals
	n = (long)ndivday + 1;
	ndiv = (long*)malloc(n*sizeof(long));
	
	//initialization
	n00 = 0;
	tc = (double)number_novalue;
	
	//loop
	do{
		tc0=tc;				
		find_sunset(n00, &n0, &n1, meteo, meteolines, horizon, horizonlines, lat, lon, ST);
		
		ndiv[0]=n0;
		for(k=1;k<=ndivday-1;k++){
			ndiv[k]=(long)(n0+k*(n1-n0)/(double)ndivday);
		}
		ndiv[ndivday]=n1;
		
		for(k=1;k<=ndivday;k++){
			tc = average_cloudiness(ndiv[k-1], ndiv[k], meteo, meteolines, lat, lon, ST, Z, sky, A);
			
			//cloudiness at night (from n00<=n<n0)
			if (k==1) {
				for (n=n00; n<n0; n++) {
					if ( (long)tc0 != number_novalue && (long)tc != number_novalue ) {
						cloudtrans[n] = tc0 + (tc-tc0) * (n-n00) / (n0-n00);
					}else {
						cloudtrans[n] = (double)number_novalue;
					}
					printf("n = %ld/%ld\n",n+1,meteolines);
				}
			}
			
			//cloudiness during the day
			for (n=ndiv[k-1]; n<ndiv[k]; n++) {
				cloudtrans[n] = tc;
				printf("n = %ld/%ld\n",n+1,meteolines);
			}
		}
		
		n00 = n1;
		
	}while (n00 < meteolines-1);
	
	n = n00;
	sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);
	alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon - ST*Pi/12. + Et)/omega);
	direction = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon - ST*Pi/12. + Et)/omega);
	if( shadows_point(horizon, horizonlines, alpha*180./Pi, direction*180./Pi, Tol_h_mount, Tol_h_flat) == 0){		
		tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, A);	
	}else {
		tc = (double)number_novalue;
	}
	cloudtrans[n] = tc;
	printf("n = %ld/%ld\n",n+1,meteolines);

}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

double find_cloudiness(long n, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double A){
	
	double tau_cloud;
	double E0, Et, Delta;
	double JD, JDbegin, JDend;
	double P, RH, T;
	double alpha, tau_atm;
	long d, m, y, h, mi;
	FILE *f;
	char *temp;
		
	//initial and final JD of the time step
	if (n == 0) {
		JDbegin = meteo[n][iJDfrom0] - 0.5 * ( meteo[n+1][iJDfrom0] - meteo[n][iJDfrom0] );
	}else {
		JDbegin = 0.5 * ( meteo[n-1][iJDfrom0] + meteo[n][iJDfrom0] );
	}
	if (n == meteolines - 1) {
		JDend = meteo[n][iJDfrom0] + 0.5 * ( meteo[n][iJDfrom0] - meteo[n-1][iJDfrom0] );
	}else {
		JDend = 0.5 * ( meteo[n][iJDfrom0] + meteo[n+1][iJDfrom0] );
	}
	
	//sun variables
	sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);
	
	//pressure [mbar]
	P=meteo[n][iPs];
	if((long)P == number_novalue || (long)P == number_absent) P=pressure(Z, 0.0, Pa0);
	
	//relative humidity [-]
	RH=meteo[n][iRh];
	if((long)RH == number_novalue || (long)RH == number_absent){
		if ( (long)meteo[n][iT] != number_absent && (long)meteo[n][iT] != number_novalue && (long)meteo[n][iTdew] != number_absent && (long)meteo[n][iTdew] != number_novalue){
			RH=RHfromTdew(meteo[n][iT], meteo[n][iTdew], Z);
		}else {
			RH=0.4;
		}
	}
	if(RH<0.01) RH=0.01;
	
	//air temperature [C]
	T=meteo[n][iT];
	if((long)T == number_novalue || (long)T == number_absent) T=0.0;
		
	//cloudiness transmissivity
	tau_cloud = cloud_transmittance(JDbegin, JDend, lat, Delta, (lon-ST*Pi/12.+Et)/omega, RH, T, P, meteo[n][iSWd],
									 meteo[n][iSWb], meteo[n][iSW], E0, sky, A);
	
	//plotting
	temp = join_strings(WORKING_DIRECTORY,filecloud);
	f = fopen(temp,"a");
	convert_JDfrom0_JDandYear(meteo[n][iJDfrom0], &JD, &y);
	convert_JDandYear_daymonthhourmin(JD, y, &d, &m, &h, &mi);
	alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon-ST*Pi/12.+Et)/omega);
	tau_atm = atm_transmittance(alpha, P, RH, T);
	fprintf(f,"%02.f/%02.f/%04.f %02.f:%02.f,%f,%f,%f,%f,%f,%f,%f\n",(float)d, (float)m, (float)y, (float)h,(float)mi, 
			alpha*180./Pi, SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon-ST*Pi/12.+Et)/omega) * 180./Pi,
			Fmax(sin(alpha), 0.05), meteo[n][iSW], Isc*E0*Fmax(sin(alpha),0.05)*tau_atm,tau_atm,tau_cloud);
	fclose(f);	
	free(temp);
	
	return tau_cloud;
		
}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

double average_cloudiness(long n0, long n1, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double A){
	
	long n;
	double tc, tc_av=0.0;
	short is_novalue=0;
		
	for(n=n0;n<n1;n++){
		tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, A);
		if( (long)tc == number_novalue){
			is_novalue = 1;
		}else {
			tc_av += tc / ((double)(n1-n0));
		}		
	}
	
	if(is_novalue==1 || n0==n1) tc_av = (double)number_novalue;
	
	return tc_av;
	
}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//nist line at which the computation begins
//n0 line at sunrise
//n1 line at sunset
//lat and lon in [rad]
void find_sunset(long nist, long *n0, long *n1, double **meteo, long meteolines, double **horizon, long horizonlines, 
				 double lat, double lon, double ST){
	
	short shad, shad0, a=0;
	long n=nist;
	double alpha, direction, E0, Et, Delta;	
	
	shad=1;//suppose it is in shadow
	*n0=-1;
		
	do{
		
		//find if it in shadow or not
		shad0=shad;		
		
		sun( meteo[n][iJDfrom0], &E0, &Et, &Delta );
		alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon - ST*Pi/12. + Et)/omega);
		direction = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon - ST*Pi/12. + Et)/omega);

		shad=shadows_point(horizon, horizonlines, alpha*180./Pi, direction*180./Pi, Tol_h_mount, Tol_h_flat);		
		
		//from shadow to non-shadow = sunrise
		if(shad0==1 && shad==0) *n0=n;
		
		//from non-shadow to shadow = sunset
		if(shad0==0 && shad==1) a=1;
		
		n++;
		
	}while(a==0 && n<meteolines);
	
	*n1=n-1;
	if(*n0==-1) *n0=(*n1);
	
}


//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
