
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "clouds.h"

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//lat and lon in [deg]
//short fill_meteo_data_with_cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat,
//									 double lon, double ST, double Z, double sky, double SWrefl_surr, long ndivday, double rotation){
//
//	double *cloudtrans;
//	long n;
//
//	//if there are radiation data, and no cloudiness
//	if ( (long)meteo[0][iSW] != number_absent || ( (long)meteo[0][iSWb] != number_absent && (long)meteo[0][iSWd] != number_absent ) ){
//
//		cloudtrans = (double*)malloc(meteolines*sizeof(double));
//		cloudiness(meteo, meteolines, horizon, horizonlines, lat*Pi/180., lon*Pi/180., ST, Z, sky, SWrefl_surr, cloudtrans, ndivday, rotation);
//
//		for (n=0; n<meteolines; n++) {
//			meteo[n][itauC] = cloudtrans[n];
//		}
//
//		free(cloudtrans);
//		return 1;
//
//	}else {
//
//		return 0;
//
//	}
//
//
//}
	
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

void cloudiness(double **meteo, long meteolines, double **horizon, long horizonlines, double lat, double lon, 
				double ST, double Z, double sky, double SWrefl_surr, double *cloudtrans, long ndivday, double rotation){
	/*!
	 * This function calculates the average cloudiness of a dataset of meteo variables recorded in a meteo station where the Short Wave radiation is measured.
	 * The day is considered to be the composed by: daylight (time where the sun is present) and night. The start of the daylight is given by moment (meteo line)
	 * in which the sun rises and the end of the daylight, corresponding with the line of the beginning of the night, corresponds to the moment (meteo line)
	 * when the sun sets.
	 * The daylight is divided into a number of intervals given by ndivday (default is 3).
	 * The function determines the meteo lines when the sun rises (n0) and when the sun sets (n1), then divides this time frame into ndivday intervals.
	 * Then, for each meteo line inside this "daylight" interval, calculates the cloud_transmittance and then makes the average inside the interval.
	 * For the cloudiness at night, linearly interpolates between the previous sunset and the sun rise.
	 * Returns a vector of average cloud_transmissivity (*cloudtrans).
	 * */

	double tc;// tau cloud
	double tc0;// tau cloud at the previous time step
	double E0;
	double Et;
	double Delta;
	double alpha; // solar elevation angle
	double direction;// solar azimuth angle
	long n00;//line at which the computation begins
	long n0;// line of the meteo file at which the sun rises
	long n1;// line of the meteo file at which the sun sets
	long k, n;// indexes
	long *ndiv;
	FILE *f;
//	char *temp;
	std::string temp;
	
	//file header
//	temp = join_strings(WORKING_DIRECTORY,filecloud);
	temp = WORKING_DIRECTORY + filecloud;
	f = fopen(temp.c_str(),"w");
	fprintf(f,"Date,SolarHeight[deg],SolarAzimuth[deg],SinSolarHeight,SWinMeasured[W/m2],SWinClearSky[W/m2],AtmTransmissivity,CloudTransmissivity\n");
	fclose(f);	
//	free(temp);
	
	//average cloudiness in ndivday daily intervals
	n = ndivday + 1;
	ndiv = (long*)malloc(n*sizeof(long));
	
	//initialization
	n00 = 0;
	tc = (double)number_novalue;
	
	//loop
	do{
		tc0=tc;				
		find_sunset(n00, &n0, &n1, meteo, meteolines, horizon, horizonlines, lat, lon, ST, rotation);
		
		ndiv[0]=n0;// line of the meteo data at which the sun rises
		for(k=1;k<=ndivday-1;k++){
			ndiv[k]=(long)(n0+k*(n1-n0)/(double)ndivday);
		}
		ndiv[ndivday]=n1;//line of the meteo data at which the sun sets
		
		for(k=1;k<=ndivday;k++){
			// find the average cloudiness during the daylight (between the line ndiv[k-1] and ndiv[k])
			tc = average_cloudiness(ndiv[k-1], ndiv[k], meteo, meteolines, lat, lon, ST, Z, sky, SWrefl_surr, rotation);
			
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
			
			//assigns the values to the vector
			for (n=ndiv[k-1]; n<ndiv[k]; n++) {
				cloudtrans[n] = tc;
				printf("n = %ld/%ld\n",n+1,meteolines);
			}
		}
		
		n00 = n1;// set the computation line at the line where the sun sets
		
	}while (n00 < meteolines-1);
	
	n = n00;
	sun(meteo[n][iJDfrom0], &E0, &Et, &Delta);
	alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega);
	direction = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega) + rotation*GTConst::Pi/180.;
	if (direction < 0) direction += 2*GTConst::Pi;
	if (direction > 2*GTConst::Pi) direction -= 2*GTConst::Pi;

	if( shadows_point(horizon, horizonlines, alpha*180./GTConst::Pi, direction*180./GTConst::Pi, GTConst::Tol_h_mount, GTConst::Tol_h_flat) == 0){
		tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, SWrefl_surr, rotation);
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

double find_cloudiness(long n, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double SWrefl_surr, double rotation){
	/*!
	 * This function calculates the cloudiness of a meteo station at a time "n" which corresponds to a line in the meteo file.
	 * According to "n", the julian day is calculated and consequently the solar parameters (elevation, azimuth).
	 * Then, the meteorological variables are read(RH, P, AirT) . Finally, the cloud transmittance is calculated through the
	 * function "cloud_transmittance" and the atmospheric transmittance through the function "atm_transmittance" that is written in the file "cloud".
	 * The function returns the cloud_transmittance.
	 *
	 * param n: line in the meteo file at which the calculation of the cloudiness is performed
	 */

	double tau_cloud;
	double E0, Et, Delta;
	double JD;// Julian Day
	double JDbegin;// julian day of begin of the calculation
	double JDend;// julian day of end of the calculation
	double P; // expected pressure at elevation Z
	double RH;
	double T;
	double alpha;
	double tau_atm;
	long d, m, y, h, mi;
	FILE *f;
//	char *temp;
	std::string temp;

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
	P=pressure(Z);
	
	//relative humidity [-]
	RH=meteo[n][iRh];
	if((long)RH != number_novalue && (long)RH != number_absent){
		RH/=100.;
	}else {
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
	tau_cloud = cloud_transmittance(JDbegin, JDend, lat, Delta, (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega, RH, T, P, meteo[n][iSWd],
									 meteo[n][iSWb], meteo[n][iSW], E0, sky, SWrefl_surr);
	
	//plotting
//	temp = join_strings(WORKING_DIRECTORY,filecloud);
	temp = WORKING_DIRECTORY + filecloud;
	f = fopen(temp.c_str(),"a");
	convert_JDfrom0_JDandYear(meteo[n][iJDfrom0], &JD, &y);
	convert_JDandYear_daymonthhourmin(JD, y, &d, &m, &h, &mi);
	alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega);
	tau_atm = atm_transmittance(alpha, P, RH, T);
	fprintf(f,"%02.f/%02.f/%04.f %02.f:%02.f,%f,%f,%f,%f,%f,%f,%f\n",(float)d, (float)m, (float)y, (float)h,(float)mi,
			alpha*180./GTConst::Pi, rotation + (SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon-ST*GTConst::Pi/12.+Et)/GTConst::omega)) * 180./GTConst::Pi,
			Fmax(sin(alpha), 0.05), meteo[n][iSW], GTConst::Isc*E0*Fmax(sin(alpha),0.05)*tau_atm,tau_atm,tau_cloud);
	fclose(f);	
//	free(temp);
	return tau_cloud;
		
}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************

double average_cloudiness(long n0, long n1, double **meteo, long meteolines, double lat, double lon, double ST, double Z, double sky, double SWrefl_surr, double rotation){
	/*! finds the average cloudiness between two times:
	 * n0: first time corresponding to a line in the meteofile
	 * n1: second time corresponding to a line in the meteofile*/
	long n;
	double tc, tc_av=0.0;
	short is_novalue=0;
		
	for(n=n0;n<n1;n++){
		tc = find_cloudiness(n, meteo, meteolines, lat, lon, ST, Z, sky, SWrefl_surr, rotation);
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

void find_sunset(long nist, long *n0, long *n1, double **meteo, long meteolines, double **horizon, long horizonlines, 
				 double lat, double lon, double ST, double rotation){
	/* !
	 * This function calculates the lines of a meteo file at which the sun rises and sets, according to the sorrounding topography given by the horizon file.
	 * param nist: line at which the computation begins
	 * param lat: latitude [rad] of the meteo station
	 * param lon: longitude [rad] of the meteo station
	 * param ST: time difference (hours) of the meteo station with respect of the UTM
	 * param rotation:
	 * param meteo: matrix of meteo data
	 * param meteolines: number of lines of the meteo dataset
	 * param horizon: horizon file of the meteo station
	 * param horizonlines: number of the horizon file lines
	 *
	 *
	 * Returns:
	 * n0: line at which the sun rises
	 * n1: line at which the sun sets
	 * */
	
	short shad, shad0, a=0;
	long n=nist;
	double alpha, direction, E0, Et, Delta;	
	
	shad=1;//suppose it is in shadow
	*n0=-1;
		
	do{
		
		//find if it in shadow or not
		shad0=shad;		
		
		sun( meteo[n][iJDfrom0], &E0, &Et, &Delta );
		alpha = SolarHeight(meteo[n][iJDfrom0], lat, Delta, (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega);
		direction = SolarAzimuth(meteo[n][iJDfrom0], lat, Delta, (lon - ST*GTConst::Pi/12. + Et)/GTConst::omega) + rotation*GTConst::Pi/180.;
		if (direction < 0) direction += 2*GTConst::Pi;
		if (direction > 2*GTConst::Pi) direction -= 2*GTConst::Pi;
		
		shad=shadows_point(horizon, horizonlines, alpha*180./GTConst::Pi, direction*180./GTConst::Pi, GTConst::Tol_h_mount, GTConst::Tol_h_flat);
		
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
