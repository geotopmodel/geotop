
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */

#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "radiation.h"
#include "meteo.h"
#include "tabs.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shortwave_radiation(long r, long c, double alpha, double direction, double E0, short shadow, double sky, 
						 double tau_cloud, double sa, double slope, double aspect, double tau_atm, double *met_data, 
						 long *met_col, double sky_st, double A, double *SWbeam, double *SWdiff, double *cosinc,
						 LONGMATRIX *nDt_shadow, LONGMATRIX *nDt_sun)

{

	double kd;
	
	if(alpha>0){ //IN THE DAYTIME

		//cosine of the incidence angle
		*cosinc=cos(slope)*sin(alpha) + sin(slope)*cos(alpha)*cos(-aspect+direction);
		//*cosinc=sin(alpha);

		//diffuse radiation
		kd=diff2glob(tau_cloud*tau_atm);
		*SWdiff=kd*Isc*E0*tau_cloud*tau_atm*sin(alpha);

		//self-shadow
		if(*cosinc<=0.0) shadow=1;
		
		//direct radiation
		if(shadow==1){
			*SWbeam=0.0;
			nDt_shadow->co[r][c]+=1;
		}else{
			*SWbeam=(1-kd)*Isc*E0*tau_cloud*tau_atm*(*cosinc);
		}
		nDt_sun->co[r][c]+=1;
			
		//printf("->kd:%f tau_cloud:%f tau_atm:%f sin(alpha):%f cosinc:%f slope:%f dir:%f aspect:%f \n",kd,tau_cloud,tau_atm,sin(alpha),*cosinc,slope*180/Pi,direction*180/Pi,aspect*180/Pi);		
		//printf("->>SWbeam:%f SWdiff:%f\n",*SWbeam,*SWdiff);
		

	}else{ //AT NIGHT

		*cosinc=0;

		*SWbeam=0.0;
		*SWdiff=0.0;

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
	 beta = 0.55*(3.912/Vis-0.01162)*(0.02472*(Vis-5)+1.132); //Vis in km (>5km)
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
	 f = fopen(files->co[ferr], "a");
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

short shadows_point(double **hor_height, double alpha, double azimuth, double tol_mount, double tol_flat)

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
	long i,buf,n=dim2(hor_height);
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

void sun(double JD, double *alpha, double *direction, double *E0, double latitude, double longitude, double standard_time){

	double Gamma, Delta, long_standard, Et, h;
	
	Gamma=2.0*Pi*JD/365.;

	//longitude of the standard time meridian 
	long_standard=standard_time*Pi/12.0;

	//correction sun-earth distance
	*E0=1.00011+0.034221*cos(Gamma)+0.00128*sin(Gamma)+0.000719*cos(2*Gamma)+0.000077*sin(2*Gamma);
	
	//Correction for sideral day (rad)	
	Et=0.000075 + 0.001868*cos(Gamma) - 0.032077*sin(Gamma) - 0.014615*cos(2*Gamma) - 0.04089*sin(2*Gamma);	/*Correction for sideral day (rad)*/
	
	//Solar Declination 
	Delta=0.006918-0.399912*cos(Gamma)+0.070257*sin(Gamma)-0.006758*cos(2*Gamma)+0.000907*sin(2*Gamma)-0.002697*cos(3*Gamma)+0.00148*sin(3*Gamma);

	//solar hour [0.0-24.0]
	h = (JD-floor(JD))*24.0 + (longitude-long_standard)/omega + Et/omega;	/*Iqbal: formula 1.4.2*/
	if(h>=24) h-=24.0;
	if(h<0) h+=24.0;
	
	/*!!NOT WORKING FOR 90 deg and -90 deg latitudes*/

	//solar height (zenith) [rad]
	*alpha=asin( sin(latitude)*sin(Delta)+cos(latitude)*cos(Delta)*cos(omega*(12-h)) );
	
	//solar azimuth (from north, clockwise) [rad]
	if(h<=12){
		if(*alpha==Pi/2.0){	//zenith
			*direction=Pi/2.0;
		}else{
			*direction=Pi - acos((sin(*alpha)*sin(latitude)-sin(Delta))/(cos(*alpha)*cos(latitude)));
		}
	}else{
		if(*alpha==Pi/2.0){ //zenith
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

	long l;
	
	/*long m;
	double res=R, z=0.0, rho, k;*/

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

void find_tau_cloud_station(long i, METEO *met, PAR *par, double alpha, double E0, double sky, double A, double *tau, double *sa){

	double P, RH, T, tau_atm, kd, kd0;
	long j;
	
	if(met->column[i-1][iPs]!=-1){
		P=met->var[i-1][met->column[i-1][iPs]];
	}else{
		P=pressure(met->st->Z->co[i], Pa0, 0.0);
	}
	
	if(met->column[i-1][iRh]!=-1){
		RH=0.01*met->var[i-1][met->column[i-1][iRh]];
	}else{
		RH=0.5;
	}
	
	if(met->column[i-1][iT]!=-1){
		T=met->var[i-1][met->column[i-1][iT]];
	}else{
		T=0.0;
	}
		
	tau_atm=atm_transmittance(alpha, P, RH, T);
	
	if(met->column[i-1][iSWb]!=-1 && met->column[i-1][iSWd]!=-1){
		if(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]>0){
			*sa=Fmax(0.01,sin(alpha));
			kd=met->var[i-1][met->column[i-1][iSWd]]/(met->var[i-1][met->column[i-1][iSWb]]+met->var[i-1][met->column[i-1][iSWd]]);
			*tau=Fmin(1.0,(met->var[i-1][met->column[i-1][iSWd]]/(Isc*E0*tau_atm*(*sa)))/(kd));
		}else{
			*tau=0.0;
		}
	
	}else if(met->column[i-1][iSW]!=-1){	
		kd=0.2;
		*sa=sin(alpha);
		if(*sa<0.05) *sa=0.05;
		if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);
					
		j=0;
		do{
			j++;
			kd0=kd;
			//SW = (1-kd(T))*Isc*T*sin + sky*Kd(T)*Isc*T*sin + (1-sky)*A*Isc*T*sin
			*tau=(met->var[i-1][met->column[i-1][iSW]]/(Isc*E0*tau_atm*(*sa)));
			if(*tau>1) *tau=1.0;				
			if(*tau<0) *tau=0.0;				
			kd=diff2glob(*tau*tau_atm);
			if(*sa<0.10) kd=(kd*(*sa-0.05)+1.0*(0.10-(*sa)))/(0.10-0.05);
		}while(fabs(kd0-kd)>0.005 && j<1000);	
				
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shadow_haiden(short point, TOPO *top, double alpha, double direction, SHORTMATRIX *shadow)
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
	int orix,oriy,rr,cc;
	long r,c;
	float sx,sy,sz,xp,yp,x,y,zray,ztopo,ct,z1,z2;
	double r2d=180.0/Pi;	//from rad to degree
	double dx=UV->U->co[1];
	double dy=UV->U->co[2];
	float MAXELEV=8848.0;
	
	initialize_shortmatrix(shadow,0); /* initialized as if it was always NOT in shadow*/
	
	if (point==1){/* PUNCTUAL (1D) simulation */
		for(r=1;r<=Nr;r++){
			for(c=1;c<=Nc;c++){
				if(top->Z0->co[r][c]!=UV->V->co[2]) shadow->co[r][c]=shadows_point(top->horizon_height[r-1][c-1], alpha*r2d, direction*r2d, 0.0, 0.0);
			}
		}

	}else{ /* DISTRIBUTED (2D) simulation */
		
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

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
