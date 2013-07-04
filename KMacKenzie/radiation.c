#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "radiation.h"
#include "shadows.h"
#include "tabs.h"
#include "meteo.09375.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, double tau_cloud, double sa, double slope, double aspect, 
	double tau_atm, double *met_data, long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc, LONGMATRIX *nDt_shadow, 
	LONGMATRIX *nDt_sun)

{

	double kd;
	
	if(alpha>0){ //DI GIORNO

		/* coseno dell'angolo di incidenza */
		*cosinc=cos(slope)*sin(alpha) + sin(slope)*cos(alpha)*cos(-aspect+direction);
		//*cosinc=sin(alpha);
		
		/* nel caso vi sia ombra propria (self shadow)*/
		if(*cosinc<=0.0) shadow=1;

		//direct and diffuse radiation
		kd=diff2glob(tau_cloud*tau_atm);
		*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
		*SWdiff=sky*kd*Isc*E0*tau_cloud*tau_atm*sin(alpha) + (1-sky)*A*Isc*E0*sin(alpha);

		//shadows
		*SWbeam*=(1-shadow);

		if(shadow==1) nDt_shadow->co[r][c]+=1;
		nDt_sun->co[r][c]+=1;
		
		/* Della Chiesa-Bertoldi   Negative Radiation check added  02.03.2010*/
		if((*SWdiff+*SWbeam)<-2){
			printf("Negative Shortwave radiation in function shortwave_radiation\n");
			printf("*SWdiff=%f, *SWbeam=%f, kd=%f, tau_cloud=%f, tau_atm=%f, *cosinc=%f\n",*SWdiff,*SWbeam,kd,tau_cloud,tau_atm,*cosinc);
			stop_execution();
		}

	}else{ /* se di notte */

		*cosinc=0;

		*SWbeam=0.0;
		*SWdiff=0.0;

	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shadow_n(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)


/*
Inputs: Z0 		matrice delle elevazioni
		curv	matrice delle curvature
		day 	giorno di osservazione
		hour	ora del giorno in questione

Outputs:alpha altezza solare in radianti
		direction: azimuth sole (da N, orario) in radianti
		shadow: matrice con le ombre (1 ombra 0 sole)
		E0:		correzione distanza Terra-Sole
*/
{
long quadrata;
double beta;
long i,j;
double r2d=180.0/Pi;	//from rad to degree

initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/

/* ###################################################################   */
if (point==1){/* PUNCTUAL (1D) simulation */
		
	for(i=1;i<=Nr;i++){
		for(j=1;j<=Nc;j++){
			if(top->Z0->co[i][j]!=UV->V->co[2]) shadow->co[i][j]=shadows_point(top->horizon_height[i-1][j-1], alpha*r2d, direction*r2d, 0.0, 0.0);
		}
	}

/* ###################################################################   */
}else{ /* DISTRIBUTED (2D) simulation */


	/**======== CALCOLO DELLE OMBRE ===============================*/
	quadrata=2*(Nc+Nr);

 	/*  Chiama Orizzonte# in geotoplib.c
		Inputs:  	U->co[1]: dim. pixel (funziona solo per pixel quadrati)
 	    			quadrata: dimensione matrice
 	    			alpha: altezza solare
 	    			Z0: matrice elevazioni
 	    			curv: matrice curvature
 	     			beta: azimuth +#Pi/4
 	     			novalue: novalue per Z0
 	    Outputs:	shadow: matrice ombre (1 ombra 0 sole) */

	if(direction>=0. && direction<=Pi/4.){
		beta=direction;
		Orizzonte1(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi/4. && direction<=Pi/2.){
		beta=(Pi/2.-direction);
		Orizzonte2(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi/2. && direction<=Pi*3./4.){
		beta=(direction-Pi/2.);
		Orizzonte3(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*3./4. && direction<=Pi){
		beta=(Pi-direction);
		Orizzonte4(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi && direction<=Pi*5./4.){
		beta=(direction-Pi);
		Orizzonte5(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*5./4. && direction<=Pi*3./2.){
		beta=(Pi*3./2.-direction);
		Orizzonte6(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*3./2. && direction<=Pi*7./4.){
		beta=(direction-Pi*3./2.);
		Orizzonte7(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);

	}else if(direction>Pi*7./4. && direction<2.*Pi){
		beta=(2.*Pi-direction);
		Orizzonte1(UV->U->co[1],quadrata,beta,alpha,top->Z0,top->curv,shadow,UV->V->co[2]);
		//error!!!
	}
} /* end of 2D simulation case*/

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//****Ratio diffuse to global radiation - Erbs et al.(1982)
double diff2glob(double a)
{
double k;
if(a<0.22){
	k=1.0-0.09*a;
}else if(a<0.80){
	k=0.9511-0.1604*a+4.388*pow(a,2.0)-16.638*pow(a,3.0)+12.336*pow(a,4.0);
}else{
	k=0.165;
}
return(k);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double atm_transmittance(double alpha, double P, double RH, double T, double A, double Vis, double Lozone)
/* Author:  Stefano Endrizzi  Year:
	* function that calculates the atmospheric transmittance
	* Input:
	* 		alpha:	solar altitude angle []
	* 		P:	air pressure
	* 		RH: air humidity [%]
	* 		T: air temperature [ûC]
	* 		A: Albedo of sorrounding terrain
	* 		Vis:	Visibility parameter (meteorological range) [Km]: gives an idea of the turbidity of the air 5<Vis<180Km
	* 		Lozone: Ozone level [cm(NTP)]: it is the height of gaseous ozone if all the ozone in a vertical column of unit area were brought
	* 				to normal temperature and surface pressure (NTP). This level changes according to latitude and seasons. Around the equator,
	* 				total ozone sverages about 0.24cm(NTP) and the amount increases with latitude, e.g. in the polar regions, it is up to 0.46cm(NTP).
	* 				It is maximum in Spring and minimum in fall (Iqbal 1983)
	* Output:
	* 		tau_atm: atmospheric transmittance
	* comment: Matteo Dall'Amico, May 2009 */
{
	double mr,//relative optical air mass
		 ma,//relative optical air mass for local condition (pressure corrected)
		 w, // pressure- and temperature-corrected precipitable water
		 beta, // Angstrom turbidity coefficient
		 g_beta,
		 alpha_wv, // absorbtance of direct irradiance by water vapor
		 alpha_g, // absorbtance of uniformly mixex gases
		 alpha_o, // ozone absorptance
		 tau_r, // trasmittance due to Rayleigh scattering
		 tau_as, // transmittance due to scattering aerosols
		 alpha_a, // absorbtance of direct irradiance by aerosols
		 tau_r_p, // trasmittance due to Rayleigh scattering when ma=1.66*P/P0
		 tau_as_p, // trasmittance due to scattering aerosols when ma=1.66*P/P0
		 tau_atm; // output

	//double p,dp;

	//transmissivity under cloudless sky (Iqbal par. 7.5)
	/* Commented 3.2.10 by Della Chiesa-Bertoldi to fix negative radiation
	mr = 1.0/(sin(alpha)+0.15*(pow((3.885+alpha*180.0/Pi),-1.253)));// eq. 5.7.2 pg. 100 Iqbal
	ma = mr*P/1013.25;// eq. 5.7.2 pg. 100 Iqbal
	w = 0.493*RH*(exp(26.23-5416.0/(T+tk)))/(T+tk); //[cm] by Leckner: see eq. 5.4.6 pg. 94 Iqbal
	beta = 0.55*(3.912/Vis-0.01162)*(0.02472*(Vis-5)+1.132); // eq. 6.6.2 pg. 119 Iqbal (5Km<Vis<180 km)
	g_beta = -0.914 + 1.909267*exp(-0.667023*beta);// eq. 7.5.8 pg. 183 Iqbal
	alpha_wv = 0.110*(pow(w*mr+6.31E-4,0.3))-0.0121;// eq. 7.5.2 pg. 182 Iqbal
	alpha_g = 0.00235*pow(126*ma+0.0129,0.26)-7.5E-4+7.5E-3*pow(ma,0.875);// eq. 7.5.3 pg. 182 Iqbal
	alpha_o = 0.045*(pow(Lozone*mr+8.34E-4,0.38))-3.1E-3;// eq. 7.5.4 pg. 182 Iqbal
	tau_r = 0.615958+0.375566*exp(-0.221185*ma);// eq. 7.5.6 pg. 183 Iqbal
	tau_as = pow(g_beta,ma);// eq. 7.5.7 pg. 183 Iqbal
	alpha_a = 0.05*tau_as;// eq. 7.5.9 pg. 183 Iqbal
	tau_r_p = 0.615958+0.375566*exp(-0.221185*1.66*P/1013.25);// eq. 7.5.14a pg. 185 Iqbal
	tau_as_p = pow(g_beta,1.66*P/1013.25);// eq. 7.5.14b pg. 185 Iqbal
	tau_atm = (1 - alpha_wv - alpha_g - alpha_o - alpha_a)*( tau_r*tau_as + 0.5*(1-tau_r) + 0.75*(1-tau_as) );// eq. 7.5.1  eq.7.5.10 and eq.7.5.11 Iqbal
	tau_atm *= (1.0 + A*(1-alpha_wv-alpha_g-alpha_o-alpha_a)*(0.5*(1-tau_r_p) + 0.25*(1-tau_as_p)));

	/*sat_vap_pressure(&p, &dp, T, P);
	tau_atm=sin(alpha)/(1.2*sin(alpha)+RH*p*(1.0+sin(alpha))*1.E-3+0.0455);*/
	
	/* Uncommented 3.2.10 by Della Chiesa to fix negative radiation  - simple Eagleson parametrization */
	tau_atm=0.47+0.47*sin(alpha);

	return(tau_atm);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min){

	double eep, fcloud;

	if(taucloud>=1){
		fcloud=0.0;
	}else{
		fcloud=pow((1.0-taucloud)/0.75,1/3.4);
		if(fcloud>1) fcloud=1.0;
	}	
			
	if(state==0){
		*eps=1.5689*pow(pvap/(T+tk),0.176)*(1.0+0.20*pow(fcloud,2.0));	//fitting of Zongo data
		*eps_min=1.5689*pow(pvap/(T+tk),0.176);
		*eps_max=*eps_min*(1.0+0.2);

	}else if(state==1){
		*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0); //Brutsaert, 1975 + Pirazzini
		//*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0+0.22*pow(fcloud,3.0)); //Brutsaert, 1975 + Pirazzini		
		*eps_min=1.24*pow((pvap/(T+tk)),0.14285714);
		//*eps_max=1.24*pow((pvap/(T+tk)),0.14285714)*1.22;
		*eps_max=0.979;
		//*eps=*eps_min*(1.0+0.44*RH-0.18*taucloud);

	}else if(state==2){
		*eps=1.08*(1.0-exp(-pow(pvap,(T+tk)/2016.0)))*(1.0+0.40*pow(fcloud,2.0));	//Satterlund, 1979
		*eps_min=1.08*(1.0-exp(-pow(pvap,(T+tk)/2016.0)));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==3){
		*eps=0.765*(1.0+0.40*pow(fcloud,2.0));	//Koenig-Langlo & Augstein, 1994
		*eps_min=0.765;
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==4){
		eep=(0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	//IDSO + HODGES + PIRAZZINI
		eep=-0.792+3.161*eep-1.573*eep*eep;
		*eps=eep*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0);	//Idso, 1981
		*eps_min=eep;
		*eps_max=0.979;

	}else if(state==5){
		*eps=(0.601+5.95*0.00001*pvap*exp(1500.0/(T+tk)))*(1.0+0.40*pow(fcloud,2.0));	//Andreas and Ackley, 1982
		*eps_min=(0.601+5.95*0.00001*pvap*exp(1500.0/(T+tk)));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==6){
		*eps=(0.23+0.484*pow(pvap/(T+tk),0.125))*(1.0+0.40*pow(fcloud,2.0));	//Konzelmann (1994)
		*eps_min=(0.23+0.484*pow(pvap/(T+tk),0.125));
		*eps_max=*eps_min*(1.0+0.4);

	}else if(state==11){
		*eps=1.24*pow((pvap/(T+tk)),0.14285714)*(1.0+0.26*fcloud); //Brutsaert, 1975+Jacobs
		*eps_min=1.24*pow((pvap/(T+tk)),0.14285714);
		*eps_max=*eps_min*(1.0+0.26);

	}else if(state==14){
		eep=(0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	//IDSO + PIRAZZINI
		*eps=eep*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0);
		*eps_min=eep;
		*eps_max=0.979;

	}else{
		t_error("Incorrect value for longwave radiation formula");
	}

	//Pirazzini et al. (2000)
	if(*eps>0.979) *eps=0.979;
	if(*eps_max>0.979) *eps_max=0.979;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double SB(double T){	//Stefan-Boltzmann law
	double R;
	R=5.67E-8*pow(T+tk,4.0);
	return(R);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double dSB_dT(double T){
	double dR_dT;
	dR_dT=4.0*5.67E-8*pow(T+tk,3.0);
	return(dR_dT);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat)

/*routine that tells you whether a point is in shadow or not, depending on the solar azimuth, elevation and horizon file at that point
 * Author: Matteo Dall'Amico, May 2008
Inputs: DOUBLEMATRIX* hor_height: matrix of horizon_height at the point
		double alpha: solar altitude (degree)
		double azimuth: solar azimuth (degree)
		double tol_mount: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
		double tol_flat: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
Output: shad: 1=the point is in shadow, 0 the point is in sun
*/
{
	double horiz_H;// horizon elevation at a defined solar azimuth
	long i,buf=0,n=dim2(hor_height); /* modified by Emanuele Cordano on 24/9/9 */
	short shad=1; //  initialized as if it was in shadow

	if(azimuth>=hor_height[n-1][0] || azimuth<hor_height[0][0]) 
		buf=1;
	/* compare the current solar azimuth with the horizon matrix */
	for (i=1; i<=n-1; i++){
		if(azimuth>=hor_height[i-1][0] && azimuth<hor_height[i][0] ) 
			buf=i+1;
		}
	horiz_H=hor_height[buf-1][1]; // horizon elevation at a particular time
	
	if(alpha<tol_flat){
		shad=1;
	}else if(alpha<horiz_H+tol_mount){
		shad=1;
	}else{
		shad=0;
	}
	
	return (shad);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void sun(double hour, double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time){

	double Gamma, Delta, long_standard, Et, h;
	
	Gamma=2.0*Pi*JD/365.;

	/* latitudine e longitudine stazione met */
	long_standard=standard_time*Pi/12.0;

	/* correzione distanza Terra-Sole */
	*E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
	/* correggo l'ora data in UTC ed ottengo l'ora reale locale*/
	Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(2*Gamma) - 0.04089*sin(2*Gamma);	/*Correction for sideral day (rad)*/
	// Solar Declination (Declinazione solare) 
	Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);

	h=hour + (longitude-long_standard)/omega + Et/omega;	/*Iqbal: formula 1.4.2*/
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;
	
	/**========= CALCOLO IL MOTO DEL SOLE ===============*/
	/*!!non funziona per latitudine=90gradi e -90gradi!!*/

	/* alpha: altezza solare in radianti */
	*alpha=asin( sin(latitude)*sin(Delta)+cos(latitude)*cos(Delta)*cos(omega*(12-h)) );
	
	/* direction: azimuth sole (da N, orario) in radianti */
	if(h<=12){
		if(*alpha==Pi/2.0){	/*sole allo zenit*/
			*direction=Pi/2.0;
		}else{
			*direction=Pi - acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}else{
		if(*alpha==Pi/2.0){ /*sole allo zenit*/
			*direction=3*Pi/2.0;
		}else{
			*direction=Pi + acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, SNOW *snow){

//	double z=0.0, res=R, rho, k;
	long l; //m;

	frac->co[1]=R;
	
	for(l=2;l<=snow->lnum->co[r][c]+1;l++){
		frac->co[l]=0.0;
	}

	/*for(l=snow->lnum->co[r][c];l>=1;l--){
		m=snow->lnum->co[r][c]-l+1;
		z+=0.001*snow->Dzl->co[l][r][c];
		rho=(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]);
		k=rho/3.0+50.0;
		//k*=1.E10;
		frac->co[m]=res-R*exp(-k*z);
		res=R*exp(-k*z);
		//printf("l:%ld res:%f R:%f k:%f %f %f %f\n",l,res,R,k,snow->w_ice->co[l][r][c],snow->w_liq->co[l][r][c],snow->Dzl->co[l][r][c]);
		if(l==1){
			frac->co[m]+=res;
			res=0.0;
		}
	}
	
	frac->co[snow->lnum->co[r][c]+1]=res;*/
	
}

	
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double A, double sky, double *tau, double *sa){
	/* Author: Stefano Endrizzi    Year:
	 * function that calculates the transmittance of the sky given measures of SW radiation
	* Input:	i:	number of meteo station under consideration
	* 			met: meteo data
	* 			par: parameters
	* 			alpha: solar altitude angle
	* 			E0: earth-sun correction
	* 			A:
	* 			sky: sky view factor on the considered station
	*
	* Output:	tau: transmittance of the sky
	* 			sa: sin(alpha)
	* comment: Matteo Dall'Amico, April 2009 */
	double P, RH, T, tau_atm, kd, kd0;

	if(met->column[i-1][iPs]!=-1){/* there is a column of air pressure in the station meteo file */
		P=met->var[i-1][met->column[i-1][iPs]];
	}else{
		P=pressure(met->st->Z->co[i], Pa0, 0.0);
	}

	if(met->column[i-1][iRh]!=-1){/* there is a column of relative humidity in the station meteo file */
		RH=0.01*met->var[i-1][met->column[i-1][iRh]];
	}else{
		RH=0.5;
	}

	if(met->column[i-1][iT]!=-1){/* there is a column of air temperature in the station meteo file */
		T=met->var[i-1][met->column[i-1][iT]];
	}else{
		T=0.0;
	}
	
	tau_atm=atm_transmittance(alpha, P, RH, T, A, par->Vis, par->Lozone);

	if(met->column[i-1][iSWb]!=-1 && met->column[i-1][iSWd]!=-1){/* there is a column of SWglobal or SWdirect in the station meteo file */
		if(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]>0){
			*sa=Fmax(0.01,sin(alpha));
			kd=met->var[i-1][met->column[i-1][iSWd]]/(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]);
			*tau=Fmin(1.0,(met->var[i-1][met->column[i-1][iSWd]]/(Isc*E0*tau_atm*(*sa)))/(sky*kd+(1-sky)*A));
		}else{
			*tau=0.0;
		}

	}else if(met->column[i-1][iSW]!=-1){/* there is a column of SWin in the station meteo file */
		kd=0.2;
		*sa=sin(alpha);
		if(*sa<0.05) *sa=0.05;
		if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);
					
		do{
			kd0=kd;
			//SW = (1-kd(T))*Isc*T*sin + sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
			*tau=(met->var[i-1][met->column[i-1][iSW]]/(Isc*E0*tau_atm*(*sa)))/((1-kd)+sky*kd+(1-sky)*A);
			if(*tau>1) *tau=1.0;				
			if(*tau<0) *tau=0.0;				
			kd=diff2glob(*tau*tau_atm);
			if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);
		}while(fabs(kd0-kd)>0.005);	
		
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

/*==================================================================================================================*/
/*==================================================================================================================*/
/*4. subroutine SHADOW_haiden*/
/*==================================================================================================================*/
/*==================================================================================================================*/
void shadow_haiden(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)
/*	Author: Thomas Haiden, Year: 16 june 2003
 * 	Function that calculates if each pixel (matrix of shadows) is in shadow or at sun given
 *  a solar elevation angle and a solar azimuth.
 *  Inputs:	top 		elevation DTM (m)
 *  		alpha:		solar elevation angle (radiants)
 *  		direction:	solar azimuth angle (from N clockwise, radiants)
 *  		point: flag indicating whether the simulation is a point (=1) or a distributed simulation
 *  Outputs:shadow: 	shadows matrix (1 shadow, 0 sun)
 *  imported by Matteo Dall'Amico on August 2009. Authorization: see email of David Whiteman on 16 june 2008
*/
{
	int orix,oriy,rr,cc;
	long r,c;
	float sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2;
	double r2d=180.0/Pi;	//from rad to degree
	double dx=UV->U->co[1];
	double dy=UV->U->co[2];
	float MAXELEV=8848.0;

	initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/

	/* ###################################################################   */
	if (point==1){/* PUNCTUAL (1D) simulation */
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->Z0->co[r][c]!=UV->V->co[2]) shadow->co[r][c]=shadows_point(top->horizon_height[r-1][c-1], alpha*r2d, direction*r2d, 0.0, 0.0);
			}
		}

	/* ###################################################################   */
	}
	else{ /* DISTRIBUTED (2D) simulation */

		/* find the sun vector components: x, y, z*/
		sx=sin(direction)*cos(alpha);
		sy=cos(direction)*sin(alpha);
		sz=sin(alpha);

		if (fabs(sx)>fabs(sy)) {
			orix=sx/fabs(sx);
			if (sy!=0) oriy=sy/fabs(sy);
			else oriy=0;
			for (r=1; r<=Nr; r++) {
				for (c=1; c<=Nc; c++) {
					rr=r;
					cc=c;
					xp=dx*cc;
					yp=dy*rr;
					zray=top->Z0->co[rr][cc];
					shadow->co[r][c]=0;// initialized at sun
					while ((shadow->co[r][c]==0)&&(zray<MAXELEV) &&(rr>1)&&(rr<=Nr-1)&&(cc>1)&&(cc<=Nc-1)) {
						ct=((cc+orix)*dx-xp)/sx;
						y=yp+ct*sy;
						if (fabs(y-dy*rr)<dy) {
							cc=cc+orix;
							xp=dx*cc;
							yp=y;
							z1=top->Z0->co[rr][cc];
							z2=top->Z0->co[rr+oriy][cc];
							ztopo=z1+(z2-z1)*(yp-dy*rr)/(oriy*dy);
						} else {
							ct=((rr+oriy)*dy-yp)/sy;
							x=xp+ct*sx;
							rr=rr+oriy;
							xp=x;
							yp=dy*rr;
							z1=top->Z0->co[rr][cc];
							z2=top->Z0->co[rr][cc+orix];
							ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx);
						}
						zray=zray+ct*sz;
						if (ztopo>zray) shadow->co[r][c]=1;
					}
				}
			}
		}
		else {
			oriy=sy/fabs(sy);
			if (sx!=0) orix=sx/fabs(sx);
			else orix=0;
			for (r=1; r<=Nr; r++) {
				for (c=1; c<=Nc; c++) {
					rr=r;
					cc=c;
					xp=dx*cc;
					yp=dy*rr;
					zray=top->Z0->co[rr][cc];
					shadow->co[r][c]=0;
					while ((shadow->co[r][c]==0)&&(zray<MAXELEV)&&(rr>1)&&(rr<=Nr-1)&&(cc>1)&&(cc<=Nc-1)) {
						ct=((rr+oriy)*dy-yp)/sy;
						x=xp+ct*sx;
						if (fabs(x-dx*cc)<dx) {
							rr=rr+oriy;
							yp=dy*rr;
							xp=x;
							z1=top->Z0->co[rr][cc];
							z2=top->Z0->co[rr][cc+orix];
							ztopo=z1+(z2-z1)*(xp-dx*cc)/(orix*dx);
						}
						else {
							ct=((cc+orix)*dx-xp)/sx;
							y=yp+ct*sy;
							cc=cc+orix;
							yp=y;
							xp=dx*cc;
							z1=top->Z0->co[rr][cc];
							z2=top->Z0->co[rr+oriy][cc];
							ztopo=z1+(z2-z1)*(yp-dy*rr)/(oriy*dy);
						}
						zray=zray+ct*sz;
						if (ztopo>zray) shadow->co[r][c]=1;
					}
				}
			}
		}
	}
}


