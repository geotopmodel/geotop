
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
#include "pedo.funct.h"
#include "constants.h"
//#include "keywords_file.h"

extern char **files;

/*
in all the following subroutine
w=theta
i=theta_ice
s=saturated water content
r=residual water content
a=alpha van genuchten
n=n van genuchten
m=m van genuchten
pmin=psi min
Ss=specific storativity
*/

/*--------------------------------------------*/
double psi_teta(double w, double i, double s, double r, double a, double n, double m, double pmin, double Ss )

{

 double psi,TETA,TETAsat,TETAmin;
 short sat;
 
 TETAsat=1.0-i/(s-r);
 TETAmin=1.0/pow((1.0+pow(a*(-pmin),n)),m);

 if(w>s-i){
	TETA=TETAsat; 
	sat=1;
 }else{
 	TETA=(w-r)/(s-r);
	sat=0;
 }
 
 if(TETA<TETAmin) TETA=TETAmin;
 
 if(TETA>1.0-1.E-6){
	psi=0.0;
 }else{
	psi=(pow((pow(TETA,-1.0/m)-1.0),1.0/n))*(-1.0/a);
 }
 if(sat==1) psi += (w-(s-i))/Ss; 

 return psi;
}




/*--------------------------------------------*/
double teta_psi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double Ss)
{
 double teta,TETA,psisat;
 short sat=0;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
 if(psi>psisat) sat=1;
  
 if(psi<pmin) psi=pmin;
 
 if(sat==0){
	if(psi>-1.E-6){
		TETA=1.0;
	}else{
		TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	}
	teta=r+TETA*(s-r);
 }else{
	teta= s-i + Ss * (psi-psisat);
 }

 return teta;
}

/*--------------------------------------------*/
double dteta_dpsi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double Ss )
//it is the derivative of teta with respect to psi [mm^-1]
{
 double dteta,psisat;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
  
 if(psi>=psisat){
	dteta=Ss;
 }else{
	dteta=(s-r)*(a*m*n)*(pow(-a*psi,n-1.0))*pow(1.0+pow(-a*psi,n),-m-1.0);
 } 

 return dteta;
}

/*--------------------------------------------*/

double T_max_dteta(double a, double n, double m){
	double T, P0, P=1.0E-3;
	long cont=0;
	do{
		P0=P;
		P=P0-d2theta(P0,n,m)/d3theta(P0,n,m);
		cont++;
	}while(fabs(P-P0)>1.E-5);
	P/=(-a);
	T=P/((1000.0*Lf)/(g*(Tfreezing+tk)));

	return(T);
}

double P_max_dteta(double a, double n, double m){
	double P0, P=1.0E-3;
	long cont=0;
	do{
		P0=P;
		P=P0-d2theta(P0,n,m)/d3theta(P0,n,m);
		cont++;
	}while(fabs(P-P0)>1.E-5);
	P/=(-a);
	return(P);
}
	
	
	
double d2theta(double P, double n, double m){
	double f, K;	
	K=n*(m+1)/(n-1);
	f=1.0/P - K*pow(P, n-1.0)/(1.0+pow(P, n));	
	return(f);
}

double d3theta(double P, double n, double m){
	double f, K;	
	K=n*(m+1)/(n-1);
	f=-1.0/(pow(P, 2.0)) - K*( (n-1.0)*pow(P, n-2.0)/(1.0+pow(P, n)) - n*pow(P, 2*n-2)/pow(1.0+pow(P, n), 2.0) );
	return(f);
}
	
	
	

/*--------------------------------------------*/
double K(double psi, double K_sat, double imp, double i, double s, double r, double a, double n, double m, double v, double pmin, double T)

{

	double TETA,psisat,K_unsat,iceratio;
	
	psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
		
	if(psi>psisat) psi=psisat;
  
	TETA=1.0/pow((1.0+pow(a*(-Fmax(pmin,psi)),n)),m);
 
	if(T>=0){
		K_unsat=K_sat*(0.000158685828*T*T+0.025263459766*T+0.731495819);
	}else{
		K_unsat=K_sat*0.73;
	}
	
	K_unsat*=(pow(TETA,v))*(pow((1-pow((1-pow(TETA,(1.0/m))),m)),2.0));
	
	iceratio=i/(s-r);
	K_unsat*=(pow(10.0,-imp*iceratio));	
							
	return K_unsat;
 
}

/*--------------------------------------------*/
double dK_dtheta(double th, double K_sat, double imp, double i, double s, double r, double a, double n, double m, double v, double pmin, double T)

{
	
	double TETA;
	double iceratio=i/(s-r);
	double dK;

	if(T>=0){
		dK = K_sat * (0.000158685828*T*T+0.025263459766*T+0.731495819);
	}else{
		dK = K_sat * 0.73;
	}
	
	dK *= (pow(10.0,-imp*iceratio));	

	TETA = (th - r)/(s - r);
	if(TETA > 0.999 - iceratio) TETA = 0.999 - iceratio;
	dK *= ( (1.-pow((1.-pow(TETA,(1./m))),m)) * pow(TETA, v-1.) * (v + 2.*TETA*pow(1.-pow(TETA, 1./m), m-1.)*pow(TETA, 1./m-1.) ) );
		
	return dK;
	
}


/*--------------------------------------------*/

double psi_saturation(double i, double s, double r, double a, double n, double m){
	
	double psisat;
	
	if(1.0-i/(s-r)>1.E-6){
		psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
	}else{
		psisat=0.0;
	}
	
	return(psisat);

}
/*------------------------------------------------------------------------------------------------------*/

double Harmonic_Mean(double D1, double D2, double K1, double K2)
	/*  calculates the harmonic mean of a property K weighted on the layer depth D
	 * It is used to calculate the property K at the interface between two consecutive layers K1,D1 and K2,D2
	 *  mean=(D1+D2)/(D1/K1+D2/K2)
	 * Author: Matteo Dall'Amico, Sept 2010 */
{
	return((D1+D2)/(D1/K1+D2/K2));
}


double Arithmetic_Mean(double D1, double D2, double K1, double K2)
{
	return((D1*K1+D2*K2)/(D1+D2));
}

/*------------------------------------------------------------------------------------------------------*/

double Mean(short a, double D1, double D2, double K1, double K2)

{

	if(a==0){
		return(Harmonic_Mean(D1,D2,K1,K2));
	}else if(a==1){
		return(Arithmetic_Mean(D1,D2,K1,K2));
	}else{
		return(0.0);
	}
}

	
/*------------------------------------------------------------------------------------------------------*/

double Psif(double T){
	
	double psi;
	
	if(T<0){
		psi=T*(1000.0*Lf)/(g*(Tfreezing+tk));
	}else{
		psi=0;
	}
		
	return(psi);
}

/******************************************************************************************************************************************/
double theta_from_psi(double psi, long l, long r, long c, SOIL *sl, double pmin){
	
	double th;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double Ss=sl->pa->co[sy][jss][l];
	
	th = teta_psi(psi, i, s, res, a, n, m, pmin, Ss);
	
	return(th);
}



/******************************************************************************************************************************************/
double psi_from_theta(double th, long l, long r, long c, SOIL *sl, double pmin){
	
	double psi;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double Ss=sl->pa->co[sy][jss][l];
	
	psi = psi_teta(th, i, s, res, a, n, m, pmin, Ss);
	
	return(psi);
}
/******************************************************************************************************************************************/
double dtheta_dpsi_from_psi(double psi, long l, long r, long c, SOIL *sl, double pmin){
	
	double dth;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double Ss=sl->pa->co[sy][jss][l];
	
	dth = dteta_dpsi(psi, i, s, res, a, n, m, pmin, Ss);
	
	return(dth);
}


/******************************************************************************************************************************************/
double k_from_psi(long jK, double psi, long l, long r, long c, SOIL *sl, double imp){

	double k;
	
	long sy=sl->type->co[r][c];

	double kmax=sl->pa->co[sy][jK][l];

	double i=sl->thice->co[l][r][c];
	double T=sl->T->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double v=sl->pa->co[sy][jv][l];
	
	k = K(psi, kmax, imp, i, s, res, a, n, m, v, PsiMin, T);
	
	
	return(k);
}

/******************************************************************************************************************************************/

double psisat_from(long l, long r, long c, SOIL *sl){

	double psi;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	
	psi = psi_saturation(i, s, res, a, n, m);
	
	return(psi);
}
/******************************************************************************************************************************************/
