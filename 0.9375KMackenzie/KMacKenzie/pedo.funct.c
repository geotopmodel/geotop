
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion KMackenzie

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 KMackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOtop is distributed in the hope that it will be useful,
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

 psi=(pow((pow(TETA,-1.0/m)-1.0),1.0/n))*(-1.0/a);

 if(sat==1) psi+=E*(w-(s-i));

 return psi;
}


/*--------------------------------------------*/
double teta_psi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E)
{
 double teta,TETA,psisat;
 short sat=0;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);
 if(psi>psisat) sat=1;

 if(psi<pmin) psi=pmin;

 if(sat==0){
	TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	teta=r+TETA*(s-r);
 }else{
	teta=s-i+(psi-psisat)/E;
 }

 return teta;
}

/*--------------------------------------------*/
double dteta_dpsi(double psi, double i, double s, double r, double a, double n, double m, double pmin, double E, double dthdp_min)
//it is the derivative of teta with respect to psi [mm^-1]
{
 double dteta,psisat;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);

 if(psi<pmin){
	dteta=dthdp_min;
 }else if(psi>psisat){
    dteta=1/E;
 }else{
	dteta=(s-r)*(a*m*n)*(pow(-a*psi,n-1.0))*pow(1.0+pow(-a*psi,n),-m-1.0);
 }

 return dteta;
}


/*--------------------------------------------*/
double K(double psi, double K_sat, double imp, double i, double s, double r, double a, double n, double m, double v, double pmin, double T)

{

 double TETA,psisat,K_unsat,iceratio;

 psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);

 if(psi>psisat) psi=psisat;
 if(psi<pmin) psi=pmin;

 TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
 iceratio=i/(r+TETA*(s-r)+i);

 if(T>=0){
	K_unsat=K_sat*(0.000158685828*T*T+0.025263459766*T+0.731495819);
 }else{
	K_unsat=K_sat*0.73;
 }
 K_unsat*=(pow(TETA,v))*(pow((1-pow((1-pow(TETA,(1.0/m))),m)),2.0));
 K_unsat*=(pow(10.0,-imp*iceratio));

 return K_unsat;
}


/*--------------------------------------------*/

double psi_saturation(double i, double s, double r, double a, double n, double m){

	double psisat;

	psisat=(pow((pow(1.0-i/(s-r),-1.0/m)-1.0),1.0/n))*(-1.0/a);

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
double Psif(double T, double psimin){

	double psi;

	if(T<0){
		psi=T*(1000.0*Lf)/(g*(Tfreezing+tk));
	}else{
		psi=0;
	}

	if(psi<psimin) psi=psimin;

	return(psi);
}

/*------------------------------------------------------------------------------------------------------*/

double Psif2(double T){

	double psi;

	if(T<0){
		psi=T*(1000.0*Lf)/(g*(Tfreezing+tk));
	}else{
		psi=0;
	}

	return(psi);
}

/*------------------------------------------------------------------------------------------------------*/

double funct_theq(double th, double W, double thi, double h, double **soil, long l, double psimin){

	double d=0.001*soil[jdz][l];
	double Cs=soil[jct][l]*(1.0-soil[jsat][l])*d;
	double A,B,C,p;
	double gg,f;

	p=(g*(Tfreezing+tk))/(1000.0*Lf);
	A=p*(Cs+W*c_ice);
	B=rho_w*d*p*(c_liq-c_ice);
	C=rho_w*d*Lf;

	gg=psi_teta(th, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);
	f=A*gg + B*th*gg + C*th - h;

	return(f);
}


/******************************************************************************************************************************************/
double dfunct_theq(double th, double W, double thi, double h, double **soil, long l, double psimin){

	double d=0.001*soil[jdz][l];
	double Cs=soil[jct][l]*(1.0-soil[jsat][l])*d;
	double A,B,C,p;
	double gg,dg,df;

	p=(g*(Tfreezing+tk))/(1000.0*Lf);
	A=p*(Cs+W*c_ice);
	B=rho_w*d*p*(c_liq-c_ice);
	C=rho_w*d*Lf;

	gg=psi_teta(th, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);
	dg=dteta_dpsi(gg, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0, 0.0);
	dg=1/dg;

	df=A*dg + B*gg + B*th*dg + C;

	return(df);
}

/******************************************************************************************************************************************/
double internal_energy_soil(double thw, double thi, double T, double D, double Csoil, double theta_sat){

	double h;

	//soil heat capacity [J m^-3 K^-1]
	//heat capacity of ice [J/(kg/K)]
	//internal energy [J/(kg*m2)]

	h=(Csoil*(1-theta_sat) + c_ice*thi*rho_w + c_liq*thw*rho_w)*T*D*0.001 + Lf*thw*rho_w*D*0.001;

	return(h);
}

/******************************************************************************************************************************************/
void from_internal_soil_energy(long r, long c, long l, double h, double *thw, double *thi, double *T, double **soil, double psimin){

	double w, th0, th1, th2, tolerance=1.0E-6;
	long cont=0, maxiter=20;
	FILE *f;

	w=( (*thw)*rho_w + (*thi)*rho_w )*soil[jdz][l]*0.001;

	//solve equation AT+BTf(kT)+Cf(kT)+D=0 - Newton Raphson
	th1=*thw;
	do{
		th0=th1;
		if(dfunct_theq(th0,w,*thi,h,soil,l,psimin)<1.0E20){
			th1=th0-funct_theq(th0,w,*thi,h,soil,l,psimin)/dfunct_theq(th0,w,*thi,h,soil,l,psimin);
		}else{
			th1=th0+10*tolerance;
		}
		//printf("CR %f %f %f %f\n",th1,th0,funct_theq(th0,w,*thi,h,soil,l,psimin),dfunct_theq(th0,w,*thi,h,soil,l,psimin));
		cont++;
	}while(cont<maxiter && fabs(th1-th0)>tolerance);

	//solve equation AT+BTf(kT)+Cf(kT)+D=0 - Bisection
	if(cont==maxiter){
		th0=soil[jres][l];
		th1=fmax(soil[jsat][l], w/soil[jdz][l]);

		if((funct_theq(th0,w,*thi,h,soil,l,psimin)>0 && funct_theq(th1,w,*thi,h,soil,l,psimin)>0)||(funct_theq(th0,w,*thi,h,soil,l,psimin)<0 && funct_theq(th1,w,*thi,h,soil,l,psimin)<0)){
			printf("\nERROR 1. IN SOIL INTERNAL ENERGY r:%ld c:%ld l:%ld\n",r,c,l);
			printf("th:%30.29f theq:%f thi:%30.29f T:%f th1:%f T1:%f\n",*thw,teta_psi(Psif2(*T), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0),*thi,*T,th1,(g*(Tfreezing+tk))/(1000.0*Lf)*psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0));
			printf("%f %f\n",funct_theq(*thw,w,*thi,h,soil,l,psimin),funct_theq(th1,w,*thi,h,soil,l,psimin));
			f=fopen(files->co[ferr]+1,"a");
			fprintf(f,"\nERROR 1. IN SOIL INTERNAL ENERGY r:%ld c:%ld l:%ld\n",r,c,l);
			fprintf(f,"th:%30.29f thi:%30.29f T:%f th1:%f T1:%f\n",*thw,*thi,*T,th1,(g*(Tfreezing+tk))/(1000.0*Lf)*psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0));
			fprintf(f,"%f %f\n",funct_theq(*thw,w,*thi,h,soil,l,psimin),funct_theq(th1,w,*thi,h,soil,l,psimin));
			fclose(f);
		}

		do{
			th2=(th0+th1)/2.0;
			if((funct_theq(th0,w,*thi,h,soil,l,psimin)>0 && funct_theq(th2,w,*thi,h,soil,l,psimin)<0)||(funct_theq(th0,w,*thi,h,soil,l,psimin)<0 && funct_theq(th2,w,*thi,h,soil,l,psimin)>0)){
				th1=th2;
			}else{
				th0=th2;
			}
		}while(fabs(th1-th0)>tolerance);

	}

	if(th1>w/soil[jdz][l]){
		*thw=w/soil[jdz][l];
		*thi=0.0;
		*T=(h-Lf*w)/( c_liq*w + soil[jct][l]*0.001*soil[jdz][l]*(1.0-soil[jsat][l]));
		//if(r==1 && c==1) printf("\nSOL1 l:%ld h:%f thw:%f thi:%f T:%f\n",l,h,*thw,*thi,*T);
	}else{
		*thw=th1;
		*thi=w/soil[jdz][l]-th1;
		*T=psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0)*(g*(Tfreezing+tk))/(1000.0*Lf);
		//if(r==1 && c==1) printf("\nSOL2 l:%ld h:%f thw:%f thi:%f T:%f funct:%f\n",l,h,*thw,*thi,*T,funct_theq(th1,w,*thi,h,soil,l,psimin));
	}

}

/******************************************************************************************************************************************/
