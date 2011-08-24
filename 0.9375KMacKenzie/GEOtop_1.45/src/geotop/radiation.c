
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

#include "constants.h"
#include "struct.geotop.h"
#include "radiation.h"
#include "meteo.h"
#include "../libraries/ascii/tabs.h"
#include "times.h"
#include "../libraries/math/util_math.h"

extern long number_novalue, number_absent;

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
//extern char *logfile;
extern long Nl, Nr, Nc;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void sun(double JDfrom0, double *E0, double *Et, double *Delta){
	
	double Gamma = 2.0*Pi*convert_JDfrom0_JD(JDfrom0)/365.25;
	
	//correction sun-earth distance
	*E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
	
	//Correction for sideral day (rad)	
	*Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(2*Gamma) - 0.04089*sin(2*Gamma);	
	
	//Solar Declination 
	*Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//angles in rad
double SolarHeight(double JD, double latitude, double Delta, double dh){
	
	double alpha;
	//solar hour [0.0-24.0]
	double h = (JD-floor(JD))*24.0 + dh;	
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;
	alpha = asin( sin(latitude)*sin(Delta) + cos(latitude)*cos(Delta)*cos(omega*(12-h)) );
	return Fmax(alpha, 0.0);
}

double SolarHeight_(double JD, double *others){
	return SolarHeight(JD, others[0], others[1], others[2]);
}

double SolarHeight__(double JD, void *others){ return SolarHeight_(JD, (double *)others); }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
//angles in rad
double SolarAzimuth(double JD, double latitude, double Delta, double dh){
	
	double alpha, direction;

	//solar hour [0.0-24.0]
	double h = (JD-floor(JD))*24.0 + dh;	
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;
	
	//solar height
	alpha=asin( sin(latitude)*sin(Delta)+cos(latitude)*cos(Delta)*cos(omega*(12-h)) );
	
	//solar azimuth
	if(h<=12){
		if(alpha==Pi/2.0){	//zenith
			direction=Pi/2.0;
		}else{
			direction=Pi - acos((sin(alpha)*sin(latitude)-sin(Delta))/(cos(alpha)*cos(latitude)));
		}
	}else{
		if(alpha==Pi/2.0){ //zenith
			direction=3*Pi/2.0;
		}else{
			direction=Pi + acos((sin(alpha)*sin(latitude)-sin(Delta))/(cos(alpha)*cos(latitude)));
		}
	}
	
	return direction;
	
}

double SolarAzimuth_(double JD, double *others){
	return SolarAzimuth(JD, others[0], others[1], others[2]);
}

double SolarAzimuth__(double JD, void *others){ return SolarAzimuth_(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double TauatmCosinc(double JD, double *others){
	
	/*double latitude = others[0];
	double Delta = others[1];
	double dh = others[2];*/
	double RH = others[3];
	double T = others[4];
	double P = others[5];
	double slope = others[6];
	double aspect = others[7];
	
	double alpha, direction;
	
	alpha = SolarHeight_(JD, others);
	if (alpha>0) {
		direction = SolarAzimuth_(JD, others);
		return atm_transmittance(Fmax(alpha,asin(0.05)),P,RH,T)*Fmax(0.0,cos(slope)*sin(alpha)+sin(slope)*cos(alpha)*cos(-aspect+direction));
	}else {
		return 0.0;
	}
}

double TauatmCosinc_(double JD, void *others){ return TauatmCosinc(JD, (double *)others); }


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double TauatmSinalpha(double JD, double *others){
	
	/*double latitude = others[0];
	double Delta = others[1];
	double dh = others[2];*/
	double RH = others[3];
	double T = others[4];
	double P = others[5];
	
	double alpha;
	
	alpha = SolarHeight_(JD, others);
	if (alpha>0) {
		return atm_transmittance(Fmax(alpha,asin(0.05)),P,RH,T) * Fmax(sin(alpha), 0.05);
	}else {
		return 0.0;
	}
}

double TauatmSinalpha_(double JD, void *others){ return TauatmSinalpha(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Cosinc(double JD, double *others){
	
	/*double latitude = others[0];
	double Delta = others[1];
	double dh = others[2];*/
	double slope = others[6];
	double aspect = others[7];
	
	double alpha, direction;
	
	alpha = SolarHeight_(JD, others);
	direction = SolarAzimuth_(JD, others);
	if (alpha>0) {
		return Fmax(0.0,cos(slope)*sin(alpha)+sin(slope)*cos(alpha)*cos(-aspect+direction));
	}else {
		return 0.0;
	}	
}

double Cosinc_(double JD, void *others){ return Cosinc(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Sinalpha(double JD, double *others){
	
	/*double latitude = others[0];
	double Delta = others[1];
	double dh = others[2];*/
	
	double alpha;
	
	alpha = SolarHeight_(JD, others);
	if (alpha>0) {
		return Fmax(sin(alpha), 0.05);
	}else {
		return 0.0;
	}
}

double Sinalpha_(double JD, void *others){ return Sinalpha(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Tauatm(double JD, double *others){
	
	/*double latitude = others[0];
	double Delta = others[1];
	double dh = others[2];*/
	double RH = others[3];
	double T = others[4];
	double P = others[5];
	
	double alpha = SolarHeight_(JD, others);
	if (alpha < asin(0.05)) alpha = asin(0.05);
	return atm_transmittance(alpha, P, RH, T);
}

double Tauatm_(double JD, void *others){ return Tauatm(JD, (double *)others); }

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shortwave_radiation(double JDbeg, double JDend, double *others, double sin_alpha, double E0, double sky, double A, 
	double tau_cloud, short shadow, double *SWb, double *SWd, double *cos_inc_bd, double *tau_atm_sin_alpha, short *SWb_yes){

	double kd, tau_atm, cos_inc, tau_atm_cos_inc;
	
	tau_atm = adaptiveSimpsons2(Tauatm_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
	//tau_atm = Tauatm( 0.5*(JDbeg+JDend), others);
	
	kd=diff2glob(tau_cloud*tau_atm);
	if(sin_alpha<0.10) kd=(kd*(sin_alpha-0.05)+1.0*(0.10-sin_alpha))/(0.10-0.05);
	
	*tau_atm_sin_alpha = adaptiveSimpsons2(TauatmSinalpha_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
	//*tau_atm_sin_alpha = tau_atm * sin_alpha;
	
	*SWd = ( Isc*E0*tau_cloud*(*tau_atm_sin_alpha) ) * ( sky*kd + (1-sky)*A );
	
	if (shadow == 1) {
		cos_inc = 0.0;
		*SWb = 0.0;
	}else {
		cos_inc = adaptiveSimpsons2(Cosinc_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
		//cos_inc = Cosinc( 0.5*(JDbeg+JDend), others);
		tau_atm_cos_inc = adaptiveSimpsons2(TauatmCosinc_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
		//tau_atm_cos_inc = tau_atm * cos_inc;
		*SWb = (1.-kd)*Isc*E0*tau_cloud*tau_atm_cos_inc;
	}
	
	*cos_inc_bd = kd*sin_alpha + (1.-kd)*cos_inc;
	
	if (sin_alpha > 1.E-5) {
		if (*SWb < 1.E-5) {
			*SWb_yes = 0;
		}else {
			*SWb_yes = 1;
		}
	}else {
		*SWb_yes = -1;
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double diff2glob(double a){

//Ratio diffuse to global radiation - Erbs et al.(1982)

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

double atm_transmittance(double X, double P, double RH, double T){
	
	//X = angle of the sun above the horizon [rad]
	//P = pressure [mbar]
	//RH = relative humidity [0-1]
	//T = air temperature [C]
	
	
	//double mr, ma, w, beta, g_beta, alpha_wv, alpha_g, alpha_o, tau_r, tau_as, alpha_a, tau_r_p, tau_as_p, tau_atm;
	//double p,dp;
	
	//FILE *f;
	
	//transmissivity under cloudless sky (Iqbal par. 7.5)
	/*mr = 1.0/(sin(alpha)+0.15*(pow((3.885+alpha*180.0/Pi),-1.253)));
	 ma = mr*P/1013.25;
	 w = 0.493*RH*(exp(26.23-5416.0/(T+tk)))/(T+tk); //cm
	 beta = 0.55*(3.912/Vis-0.01162)*(0.02472*(Vis-5)+1.1452); //Vis in km (>5km)
	 g_beta = -0.914 + 1.909267*exp(-0.667023*beta);
	 alpha_wv = 0.110*(pow(w*mr+6.31E-4,0.3))-0.0121;
	 alpha_g = 0.00235*pow(126*ma+0.0129,0.26)-7.5E-4+7.5E-3*pow(ma,0.875);
	 alpha_o = 0.045*(pow(Lozone*mr+8.34E-4,0.38))-3.1E-3;
	 tau_r = 0.615958+0.375566*exp(-0.221185*ma);
	 tau_as = pow(g_beta,ma);
	 alpha_a = 0.05*tau_as;
	 tau_r_p = 0.615958+0.375566*exp(-0.221185*1.66*P/1013.25);
	 tau_as_p = pow(g_beta,1.66*P/1013.25);
	 tau_atm = (1 - alpha_wv - alpha_g - alpha_o - alpha_a)*( tau_r*tau_as + 0.5*(1-tau_r) + 0.75*(1-tau_as) );
	 tau_atm *= (1.0 + A*(1-alpha_wv-alpha_g-alpha_o-alpha_a)*(0.5*(1-tau_r_p) + 0.25*(1-tau_as_p)));
	 
	 if (tau_atm < 0 || tau_atm > 1){
	 f = fopen(logfile, "a");
	 fprintf(f, "tau_atm:%f alpha_wv:%f alpha_g:%f alpha_o:%f tau_r:%f tau_as:%f alpha_a:%f tau_r_p:%f tau_as_p:%f\n",tau_atm,alpha_wv,alpha_g,alpha_o,tau_r,tau_as,alpha_a,tau_r_p,tau_as_p);
	 fprintf(f, "For these set of parameter the Iqbal parameterization of solar irradiance returns values out of range [0,1]. Check the Iqbal parameterization, or please refer the problem to the GEOtop team\n");
	 
	 printf("tau_atm:%f alpha_wv:%f alpha_g:%f alpha_o:%f tau_r:%f tau_as:%f alpha_a:%f tau_r_p:%f tau_as_p:%f\n",tau_atm,alpha_wv,alpha_g,alpha_o,tau_r,tau_as,alpha_a,tau_r_p,tau_as_p);
	 t_error("For these set of parameter the Iqbal parameterization of solar irradiance returns values out of range [0,1]. Check the Iqbal parameterization, or please refer the problem to the GEOtop team");	
	 }*/
	
	//other formulations
	
	//sat_vap_pressure(&p, &dp, T, P);
	//tau_atm=sin(alpha)/(1.2*sin(alpha)+RH*p*(1.0+sin(alpha))*1.E-3+0.0455);
	
	//tau_atm=0.47+0.47*sin(alpha);
	
	//from Mayers and Dale, Predicting Daily Insolation with Hourly Cloud Height and Coverage, 1983, pg 537, Journal of Climate and Applied Meteorology
	double tau_sa;//Reyleigh scattering and gas absorption transmittance
	double tau_w;//transmittance due to water vapor
	double tau_a;//transmittance due to aerosol
	double tau_atm;//global clear sky atmospheric transmittance
	double m;//optical air mass at sea level pressure
	double w;//precipitable water [cm]
	
	m = 35. * pow( 1224.*pow(sin(X), 2.) + 1. , -0.5 );
	tau_sa = 1.021 - 0.084 * pow(m*(0.000949*P + 0.051), 0.5);
	w = 0.493*RH*(exp(26.23-5416.0/(T+tk)))/(T+tk); 
	tau_w = 1. - 0.077*pow(w*m, 0.3);
	tau_a = pow(0.935, m);
	
	tau_atm = tau_sa*tau_w*tau_a;
		
	return(tau_atm);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void longwave_radiation(short state, double pvap, double RH, double T, double taucloud, double *eps, double *eps_max, double *eps_min){

	double taucloud_overcast=0.29;//after Kimball(1928)
				
	if(state==1){
		*eps_min = 1.24*pow((pvap/(T+tk)),1./7.); //Brutsaert, 1975

	}else if(state==2){
		*eps_min = 1.08*(1.0-exp(-pow(pvap,(T+tk)/2016.0)));	//Satterlund, 1979
		
	}else if(state==3){
		*eps_min = (0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	//Idso(1981)
		
	}else if(state==4){
		*eps_min = (0.7+5.95*0.00001*pvap*exp(1500/(T+tk)));	
		*eps_min = -0.792 + 3.161*(*eps_min) - 1.573*(*eps_min)*(*eps_min);	//IDSO + HODGES		

	}else if(state==5){
		*eps_min = 0.765;	//Koenig-Langlo & Augstein, 1994

	}else if(state==6){
		*eps_min = (0.601+5.95*0.00001*pvap*exp(1500.0/(T+tk)));//Andreas and Ackley, 1982

	}else if(state==7){
		*eps_min = (0.23+0.484*pow(pvap/(T+tk),0.125));	//Konzelmann (1994)
		
	}else if(state==8){
		*eps_min = (1.-(1.+46.5*pvap/(T+tk))*exp(-pow(1.2+3.*46.5*pvap/(T+tk) , 0.5)));//Prata 1996
		
	}else if(state==9){
		*eps_min = ( 59.38 + 113.7*pow( (T+tk)/273.16 , 6. ) + 96.96*pow((465.*pvap/(T+tk))/25., 0.5) ) / (5.67E-8*pow(T+tk, 4.));//Dilley 1998

	}else{
		t_error("Incorrect value for longwave radiation formula");
	}
	
	*eps = (*eps_min) * taucloud + 1.0 * (1.-taucloud);
	*eps_max = (*eps_min) * taucloud_overcast + 1.0 * (1.-taucloud);
		
	/*double fcloud;
	fcloud=pow((1.0-taucloud)/0.75,1/3.4);
	*eps=(*eps_min)*(1.0-pow(fcloud,6.0))+0.979*pow(fcloud,4.0); Pirazzini*/


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

void rad_snow_absorption(long r, long c, DOUBLEVECTOR *frac, double R, STATEVAR_3D *snow){

	long l, m;
	double res=R, z=0.0, rho, k;

	initialize_doublevector(frac, 0.);
	
	//in case of snow
	if( snow->lnum->co[r][c] > 1){
		
		for(l=snow->lnum->co[r][c];l>=1;l--){
			m=snow->lnum->co[r][c]-l+1;
			z+=0.001*snow->Dzl->co[l][r][c];
			rho=(snow->w_ice->co[l][r][c]+snow->w_liq->co[l][r][c])/(0.001*snow->Dzl->co[l][r][c]);
			k=rho/3.0+50.0;
			frac->co[m]=res-R*exp(-k*z);
			res=R*exp(-k*z);
		}
	
		frac->co[snow->lnum->co[r][c]+1]=res;
		
	}else {
		
		frac->co[0] = res;
		
	}
}
	
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double cloud_transmittance(double JDbeg, double JDend, double lat, double Delta, double dh, double RH, double T,
						   double P, double SWd, double SWb, double SW, double E0, double sky, double A){
	
	double *others;
	double tau_atm, tau_atm_sin_alpha, sin_alpha, kd, kd0;
	double tau = (double)number_novalue;
	long j;
	
	others = (double*)malloc(6*sizeof(double));
	others[0] = lat;
	others[1] = Delta;
	others[2] = dh;
	others[3] = RH;
	others[4] = T;
	others[5] = P;
	
	tau_atm_sin_alpha = adaptiveSimpsons2(TauatmSinalpha_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
	//tau_atm = Tauatm( 0.5*(JDbeg+JDend), others);
	//sin_alpha = Sinalpha( 0.5*(JDbeg+JDend), others);
	//tau_atm_sin_alpha = tau_atm*sin_alpha;
	
	if (tau_atm_sin_alpha > 0) {
		if((long)SWd!=number_absent && (long)SWd!=number_novalue && (long)SWb!=number_absent && (long)SWb!=number_novalue){
			if( SWb+SWd > 0 && SWd > 0){
				kd = SWd / (SWb+SWd);
				tau = SWd / ( Isc*E0*tau_atm_sin_alpha ) / ( sky*kd + (1-sky)*A );
			}
			
		}else if((long)SW!=number_absent && (long)SW!=number_novalue){
			
			tau_atm = adaptiveSimpsons2(Tauatm_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
			sin_alpha = adaptiveSimpsons2(Sinalpha_, others, JDbeg, JDend, 1.E-4, 10) / (JDend - JDbeg);
			
			kd=0.2;
			if(sin_alpha<0.10) kd=(kd*(sin_alpha-0.05)+1.0*(0.10-sin_alpha))/(0.10-0.05);
			
			j=0;
			do{
				
				j++;
				kd0=kd;
				//SW = (1-kd(T))*Isc*T*sin + sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
				//T=Ta*Tc
				//SW = Tc * (Isc*Ta*sin) * ( (1-kd) + sky*kd + (1-sky)*A )
				//Tc = SW / ( (Isc*Ta*sin) * ( (1-kd) + sky*kd + (1-sky)*A ) )
				tau = SW / ( Isc*E0*tau_atm_sin_alpha * ( (1-kd) + sky*kd + (1-sky)*A ) );
				if(tau > 1) tau = 1.0;				
				if(tau < 0) tau = 0.0;				
				kd = diff2glob(tau * tau_atm);
				if(sin_alpha<0.10) kd=(kd*(sin_alpha-0.05)+1.0*(0.10-sin_alpha))/(0.10-0.05);
				
			}while(fabs(kd0-kd)>0.005 && j<1000);
		}
	}
		
	if( (long)tau != number_novalue){
		if(tau<min_tau_cloud) tau=min_tau_cloud;
	}
	
	free(others);
	
	return tau;
}
	
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_tau_cloud_station(double JDbeg, double JDend, long i, METEO *met, double Delta, double E0, double Et, 
							  double ST, double A){

	double P, RH, T;
	
	//pressure [mbar]
	P=met->var[i-1][iPs];
	if((long)P == number_novalue || (long)P == number_absent) P=pressure(met->st->Z->co[i], 0.0, Pa0);
	
	//relative humidity 
	RH=met->var[i-1][iRh];
	if((long)RH == number_novalue || (long)RH == number_absent){
		if ( (long)met->var[i-1][iT] != number_absent && (long)met->var[i-1][iT] != number_novalue && (long)met->var[i-1][iTdew] != number_absent && (long)met->var[i-1][iTdew] != number_novalue){
			RH=RHfromTdew(met->var[i-1][iT], met->var[i-1][iTdew], met->st->Z->co[i]);
		}else {
			RH=0.4;
		}
	}
	if(RH<0.01) RH=0.01;
		
	T=met->var[i-1][iT];
	if((long)T == number_novalue || (long)T == number_absent) T=0.0;
			
	return cloud_transmittance(JDbeg, JDend, met->st->lat->co[i]*Pi/180., Delta, (met->st->lon->co[i]*Pi/180. - ST*Pi/12. + Et)/omega, RH,
							   T, P, met->var[i-1][iSWd], met->var[i-1][iSWb], met->var[i-1][iSW], E0, met->st->sky->co[i], A);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short shadows_point(double **hor_height, long hor_lines, double alpha, double azimuth, double tol_mount, double tol_flat)

/*routine that tells you whether a point is in shadow or not, depending on the solar azimuth, elevation and horizon file at that point
 * Author: Matteo Dall'Amico, May 2010
 Inputs: DOUBLEMATRIX* hor_height: matrix of horizon_height at the point
 double alpha: solar altitude (degree)
 double azimuth: solar azimuth (degree)
 double tol_mount: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
 double tol_flat: tolerance over a mountaneaus horizon to have a reliable cloud datum (degree)
 Output: shad: 1=the point is in shadow, 0 the point is in sun
 */
{
	double horiz_H;// horizon elevation at a defined solar azimuth
	long i,buf=1,n=hor_lines;
	short shad=1; //  initialized as if it was in shadow
	
	/* compare the current solar azimuth with the horizon matrix */
	if(azimuth>=hor_height[n-1][0] || azimuth<hor_height[0][0]){
		buf=1;
	}else{
		for (i=1; i<=n-1; i++){
			if(azimuth>=hor_height[i-1][0] && azimuth<hor_height[i][0]) buf=i+1;
		}
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

void shadow_haiden(TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)
/*	Author: Thomas Haiden, Year: 16 june 2003
 * 	Function that calculates if each pixel (matrix of shadows) is in shadow or at sun given
 *  a solar elevation angle and a solar azimuth.
 *  Inputs:	top 		elevation DTM (m)
 *  		alpha:		solar elevation angle (radiants)
 *  		direction:	solar azimuth angle (from N clockwise, radiants)
 *  		point: flag indicating whether the simulation is a point (=1) or a distributed simulation
 *  Outputs:shadow: 	shadows matrix (1 shadow, 0 sun)
 *  imported by Matteo Dall'Amico on August 2009. Authorization: see email of David Whiteman on 16 june 2010
 */
{
	long orix,oriy,rr,cc;
	long r,c;
	float sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2;
	double dx=UV->U->co[1];
	double dy=UV->U->co[2];
	float MAXELEV=8848.0;
	
	initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/
			
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double find_albedo(double dry_albedo, double sat_albedo, double wat_content, double residual_wc, double saturated_wc){

	return (dry_albedo + (sat_albedo-dry_albedo) * (wat_content - residual_wc) / (saturated_wc - residual_wc) );

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/


void find_actual_cloudiness(double *tau_cloud, double *tau_cloud_av, short *tau_cloud_yes, short *tau_cloud_av_yes, 
							METEO *met, double JDb, double JDe, double Delta, double E0, double Et, double ST, double A){
	
	//SWdata = flag -> 0= no radiation data measured, 1=beam and diff measured, 2=global measured
	short SWdata;
	double tc;
	
	if((long)met->var[met->nstsrad-1][iSWb]!=number_absent && (long)met->var[met->nstsrad-1][iSWd]!=number_absent){
		if((long)met->var[met->nstsrad-1][iSWb]!=number_novalue && (long)met->var[met->nstsrad-1][iSWd]!=number_novalue){
			SWdata=2;
		}else{
			SWdata=0;
		}
	}else if((long)met->var[met->nstsrad-1][iSW]!=number_absent){
		if((long)met->var[met->nstsrad-1][iSW]!=number_novalue){
			SWdata=1;
		}else{
			SWdata=0;
		}
	}else{
		SWdata=0;
	}
	
	if(SWdata>0){
		tc = find_tau_cloud_station(JDb, JDe, met->nstsrad, met, Delta, E0, Et, ST, A);
		if ( (long)tc != number_novalue){
			*tau_cloud_yes = 1;
			*tau_cloud = tc;
		}else {
			*tau_cloud_yes = 0;
		}		
	}else{
		*tau_cloud_yes = 0;
	}
	
	if( (long)met->var[met->nstcloud-1][iC]!=number_absent ){
		
		tc = met->var[met->nstcloud-1][iC];
		
		if((long)tc!=number_novalue){
			*tau_cloud_av_yes = 1;
			tc = 1. - 0.71*tc;//from fraction of sky covered by clouds to cloud transmissivity after Kimball (1928)
			if(tc > 1) tc = 1.;
			if(tc < 0) tc = 0.;	
			*tau_cloud_av = tc;
		}else{
			*tau_cloud_av_yes = 0;
		}
		
	}else if( (long)met->var[met->nstcloud-1][itauC]!=number_absent ){
		
		tc = met->var[met->nstcloud-1][itauC];
		
		if((long)tc!=number_novalue){
			*tau_cloud_av_yes = 1;
			if(tc > 1) tc = 1.;
			if(tc < 0) tc = 0.;	
			*tau_cloud_av = tc;
		}else{
			*tau_cloud_av_yes = 0;
		}
		
	}else{
		
		*tau_cloud_av_yes = 0;
		
	}
}

//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************
//*****************************************************************************************************************



