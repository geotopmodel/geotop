
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


#include "keywords_file.h"
#include "struct.geotop.09375.h"
#include "pedo.funct.h"
#include "constant.h"

extern STRINGBIN *files;

/*in all the following subroutine
w=theta
i=theta_ice*/

/*--------------------------------------------*/
double psi_teta(double w, double i, double s, double r, double a, double n, double m, double pmin, double E)

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
 if(sat==1) psi+=E*(w-(s-i));

 return psi;
}




double psi_teta2(double w, double i, double s, double r, double a, double n, double m, double pmin, double E)

{

 double psi,TETA,TETAsat,TETAmin;
 short sat;

 TETAsat=1.0-i/(s-r);
 TETAmin=1.0/pow((1.0+pow(a*(-pmin),n)),m);

 printf("tetasat:%e i:%e s:%e r:%e tetamin:%f\n",TETAsat,i,s,r,TETAmin);

 if(w>s-i){
	TETA=TETAsat;
	sat=1;
	printf("1: theta:%f sat:%d\n",TETA,sat);
 }else{
 	TETA=(w-r)/(s-r);
	sat=0;
	printf("2: theta:%f sat:%d\n",TETA,sat);
 }

 if(TETA<TETAmin) TETA=TETAmin;

 if(TETA>=1.0){
	printf("a\n");
	psi=0.0;
 }else{
	printf("b\n");
	psi=(pow((pow(TETA,-1.0/m)-1.0),1.0/n))*(-1.0/a);
 }

 printf("psi:%f %f %f\n",psi,pow(TETA,-1.0/m), pow((pow(TETA,-1.0/m)-1.0),1.0/n)  );

 if(sat==1) psi+=E*(w-(s-i));

 printf("E:%e i:%e s-i:%e w:%e psi:%e\n",E,i,s-i,w,psi);

 return psi;
}


/*--------------------------------------------*/
double teta_psi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E)
{
 double teta,TETA,psisat;
 short sat=0;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);// was like this
 //psisat=Fmax(pmin, (pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a));// non va bene 9/11/09
 if(psi>psisat) sat=1;

 if(psi<pmin) psi=pmin;

 if(sat==0){
	if(psi>-1.E-6){
		TETA=1.0;
		//if(fabs(psi-pmin)<1.0E-6) {printf("\nteta=%f, TETA=%F, i=%f, s=%f, r=%f",teta, TETA, i, s, r); stop_execution();}
	}else{
		TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	}
	teta=r+TETA*(s-r);
 }else{
	teta=s-i+(psi-psisat)/E;
	//if(fabs(psi-pmin)<1) {printf("\nsat!=0: teta=%f, psi=%f, pmin=%f, psisat=%f, TETA=%F, ice=%f, theta_s=%f, theta_r=%f",teta, psi, pmin, psisat,TETA, i, s, r); stop_execution();}
 }

 return teta;
}

double teta_psi2(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E)
{
 double teta,TETA,psisat;
 short sat=0;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
 if(psi>psisat) sat=1;

  printf("....psisat:%e %e %e %e %e\n",psisat,1.0-i/(s-r),i,s,r);


 if(psi<pmin) psi=pmin;

 if(sat==0){
	if(psi>-1.E-6){
		TETA=1.0;
	}else{
		TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	}
	teta=r+TETA*(s-r);
 }else{
	teta=s-i+(psi-psisat)/E;
	printf("...%f %f %f %f %f %e\n",teta,s,i,psi,psisat,E);
 }

 return teta;
}

/*--------------------------------------------*/
double dteta_dpsi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E)
//it is the derivative of teta with respect to psi [mm^-1]
{
 double dteta,psisat;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);

 if(psi>=psisat){
	dteta=1/E;
 }else{
	dteta=(s-r)*(a*m*n)*(pow(-a*psi,n-1.0))*pow(1.0+pow(-a*psi,n),-m-1.0);
 }

 return dteta;
}


double T_max_dteta(double a, double n, double m){
	double T, P0, P=1.0E-3;
	long cont=0;
	do{
		P0=P;
		P=P0-d2theta(P0,n,m)/d3theta(P0,n,m);
		cont++;
//		printf("P:%f P0:%f cont:%ld\n",P,P0,cont);
	}while(fabs(P-P0)>1.E-5);
	P/=(-a);
	T=P/((1000.0*Lf)/(g*(Tfreezing+tk)));
//	printf("P:%f T:%f\n",P,T);
//	stop_execution();
	return(T);
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

 double TETA,psisat,K_unsat,iceratio;//TETATOT;

 if(psi<pmin){
	K_unsat=0.0;
 }else{

	psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);

	if(psi>psisat) psi=psisat;

	TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	iceratio=i/(r+TETA*(s-r)+i);

	if(T>=0){
		K_unsat=K_sat*(0.000158685828*T*T+0.025263459766*T+0.731495819);
	}else{
		K_unsat=K_sat*0.73;
	}

	K_unsat*=(pow(TETA,v))*(pow((1-pow((1-pow(TETA,(1.0/m))),m)),2.0));
	K_unsat*=(pow(10.0,-imp*iceratio));

	//TETATOT=i/(s-r)+TETA;
	//if(TETATOT>0.95 && i/(s+r)>0) K_unsat*=(pow(10.0,-7.0*((TETATOT-0.95)/0.05)));
	//if(i>0) K_unsat=0.0;

 }
 return K_unsat;

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
	 * Author: Matteo Dall'Amico, Sept 2008 */
{
	return((D1+D2)/(D1/K1+D2/K2));
	//return((D1*K1+D2*K2)/(D1+D2));
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

/*------------------------------------------------------------------------------------------------------*/

double funct_T(double T, double W, double h, double **soil, long l, double psimin){

	double d=0.001*soil[jdz][l];
	double Cs=soil[jct][l]*(1.0-soil[jsat][l])*d;
	double A,B,C,p=(1000.0*Lf)/(g*(Tfreezing+tk));
	double gg,f;

	//(c_ice*W_ice + c_liq*W_liq + Cs)*(T-Tf) + Lf*W_liq - h = 0
	//(c_ice*W - c_ice*W_liq + c_liq*W_liq + Cs)*(T-Tf) + Lf*W_liq - h = 0
	//(c_ice*W + W_liq*(c_liq-c_ice) + Cs)*T + Lf*W_liq - h = 0
	//
	//Wliq = D*rho_w*teta_psi(psi(T)) = D*rho_w*teta_psi(T*k)
	//k = (1000.0*Lf)/(g*(Tfreezing+tk))
	//(c_ice*W + Cs)*T + (c_liq-c_ice)*D*rho_w*teta_psi(T*k)*T + Lf*D*rho_w*teta_psi(T*k) - h = 0
	//A*T + B*teta_psi(T*k)*T + C*teta_psi(T*k) - h = 0 = funct_T
	//
	//where
	//
	//A=c_ice*W + Cs
	//B=(c_liq-c_ice)*D*rho_w
	//C=Lf*D*rho_w

	A=Cs+W*c_ice;
	B=rho_w*d*(c_liq-c_ice);
	C=rho_w*d*Lf;

	gg=teta_psi(Fmin(p*T,0.0), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);

	//printf("%ld gg:%f psi:%f %f T:%f %f\n",l,gg,p*T,Fmin(p*T,0.0),T,soil[jsat][l]);

	f=A*T + B*gg*T + C*gg - h;

	return(f);
}


/******************************************************************************************************************************************/
double dfunct_T(double T, double W, double h, double **soil, long l, double psimin){

	double d=0.001*soil[jdz][l];
	double Cs=soil[jct][l]*(1.0-soil[jsat][l])*d;
	double A,B,C,p=(1000.0*Lf)/(g*(Tfreezing+tk));
	double gg,dg,df;

	//A*T + B*teta_psi(T*k)*T + C*teta_psi(T*k) - h = 0 = funct_T
	//first derivative
	//A + B*(teta_psi(T*k) + T*k*dteta_psi(T*k)) + C*k*dteta_psi(T*k) = dfunct_T
	//

	A=Cs+W*c_ice;
	B=rho_w*d*(c_liq-c_ice);
	C=rho_w*d*Lf;

	gg=teta_psi(Fmin(p*T,0.0), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);
	dg=dteta_dpsi(Fmin(p*T,0.0), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);

	df=A + B*(gg + T*p*dg) + C*p*dg;

	return(df);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
double theta_from_psi(double psi, long l, long r, long c, SOIL *sl, double Esoil){
	
	double th;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	
	th = teta_psi(psi, i, s, res, a, n, m, PSImin, Esoil);
	
	return(th);
}
/******************************************************************************************************************************************/
double psi_from_theta(double th, long l, long r, long c, SOIL *sl, double Esoil){
	
	double psi;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	
	psi = psi_teta(th, i, s, res, a, n, m, PSImin, Esoil);
	
	return(psi);
}
/******************************************************************************************************************************************/
double dtheta_dpsi_from_psi(double psi, long l, long r, long c, SOIL *sl, double Esoil){
	
	double dth;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	
	dth = dteta_dpsi(psi, i, s, res, a, n, m, PSImin, Esoil);
	
	return(dth);
}
/******************************************************************************************************************************************/
double k_from_psi(long jK, double psi, long l, long r, long c, SOIL *sl, double imp){

	double k;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double T=sl->T->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double kmax=sl->pa->co[sy][jK][l];
	double v=sl->pa->co[sy][jv][l];
	
	k = K(psi, kmax, imp, i, s, res, a, n, m, v, PSImin, T);
	
	return(k);
}
double k_from_psi2(long jK, double psi, long l, long r, long c, SOIL *sl, double imp){

	double k;
	
	long sy=sl->type->co[r][c];
	
	double i=sl->thice->co[l][r][c];
	double T=sl->T->co[l][r][c];
	double s=sl->pa->co[sy][jsat][l];
	double res=sl->pa->co[sy][jres][l];
	double a=sl->pa->co[sy][ja][l];
	double n=sl->pa->co[sy][jns][l];
	double m=1.-1./n;
	double kmax=sl->pa->co[sy][jK][l];
	double v=sl->pa->co[sy][jv][l];
	
	printf("%f %f %f %f %f %f %f %f %f %f %f %f\n",psi, kmax, imp, i, s, res, a, n, m, v, PSImin, T);
	
	k = K(psi, kmax, imp, i, s, res, a, n, m, v, PSImin, T);
	
	printf("k:%f\n",k);
	
	stop_execution();
	
	
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
