
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
 
 printf("tetasat:%f tetamin:%f\n",TETAsat,TETAmin);

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
 
 if(TETA>1.0-1.E-6){
	printf("a\n");
	psi=0.0;
 }else{
	printf("b\n");
	psi=(pow((pow(TETA,-1.0/m)-1.0),1.0/n))*(-1.0/a);
 }
 
 printf("psi:%f %f %f\n",psi,pow(TETA,-1.0/m), pow((pow(TETA,-1.0/m)-1.0),1.0/n)  );

 if(sat==1) psi+=E*(w-(s-i)); 

 printf("E:%f s-i:%f w:%f psi:%f\n",E,s-i,w,psi);

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
	if(psi>-1.E-6){
		TETA=1.0;
	}else{
		TETA=1.0/pow((1.0+pow(a*(-psi),n)),m);
	}
	teta=r+TETA*(s-r);
 }else{
	teta=s-i+(psi-psisat)/E;
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

 double TETA,psisat,K_unsat,iceratio;

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
	
	double thtot, T1, T0, tolerance=1.0E-8;
	long cont=0, maxiter=100;
	double n,psim,Tstar;
	double p=(1000.0*Lf)/(g*(Tfreezing+tk));
	
	thtot = (*thw) + (*thi);
										
	//solve equation AT+BTf(kT)+Cf(kT)+D=0 - Newton Raphson
	T1 = *T;
	do{
		T0 = T1;
		n = 1.0;
		cont++;
		
		do{			
			T1 = T0 - n*funct_T(T0, thtot*soil[jdz][l], h, soil, l, psimin)/dfunct_T(T0, thtot*soil[jdz][l], h, soil, l, psimin);
			//printf("..T0:%f f:%f df:%f\n",T0,funct_T(T0, thtot*soil[jdz][l], h, soil, l, psimin),dfunct_T(T0, thtot*soil[jdz][l], h, soil, l, psimin));
			n/=2.0;
		}while(fabs(funct_T(T0, thtot*soil[jdz][l], h, soil, l, psimin)) < fabs(funct_T(T1, thtot*soil[jdz][l], h, soil, l, psimin)) );

	}while(cont<maxiter && fabs(T1-T0)>tolerance);
	
	psim=psi_teta(thtot, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);
	Tstar=Fmin(psim/p, 0.0);
	
	if(T1 <= Tstar){
		*T = T1;
		*thw = teta_psi(p*(*T), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0);
		*thi = thtot - (*thw);
	}else{
		*thw = thtot;
		*thi = 0.0;
		*T = (h-Lf*(*thw)*soil[jdz][l]) / ( soil[jct][l]*(1.0-soil[jsat][l])*0.001*soil[jdz][l] + c_liq*(*thw)*soil[jdz][l] );
	}
	
	if(*T!=(*T) || (*thw)!=(*thw) || (*thi)!=(*thi)) printf("No value in soil energy balance T:%f thw:%f thi:%f l:%ld r:%ld c:%ld\n",*T,*thw,*thi,l,r,c);
	if(fabs(internal_energy_soil(*thw, *thi, *T, soil[jdz][l], soil[jct][l], soil[jsat][l]) - h ) > 0.1){
		printf("not converging soil energy balance: h:%f T:%f thw:%f thi:%f hnew:%f l:%ld r:%ld c:%ld\n",h,*T,*thw,*thi,
			internal_energy_soil(*thw, *thi, *T, soil[jdz][l], soil[jct][l], soil[jsat][l]) ,l,r,c);
	}
		
}

/******************************************************************************************************************************************/
/*void from_internal_soil_energy2(long r, long c, long l, double h, double *thw, double *thi, double *T, double **soil, double psimin){
	
	double w, th0, th1, th2, tolerance=1.0E-6;
	long cont=0, maxiter=100;
	double n;
	FILE *f;
				
	w=( (*thw)*rho_w + (*thi)*rho_w )*soil[jdz][l]*0.001;
					
	//solve equation AT+BTf(kT)+Cf(kT)+D=0 - Newton Raphson
	th1=*thw;
	do{
		th0=th1;
		n=1.0;
		cont++;
		do{
			if(th0==soil[jsat][l]) th0=soil[jsat][l]-0.001;
			if(dfunct_theq(th0,w,*thi,h,soil,l,psimin)!=0) th1=th0-n*funct_theq(th0,w,*thi,h,soil,l,psimin)/dfunct_theq(th0,w,*thi,h,soil,l,psimin);
			if(r==19 && c==19) printf("l:%ld h:%f th0:%f th1:%f thi:%f T:%f n:%f funct_theq(th0,w,*thi,h,soil,l,psimin):%f funct_theq(th1,w,*thi,h,soil,l,psimin):%f dfunct_theq(th0,w,*thi,h,soil,l,psimin):%f\n",
				l,h,th0,th1,*thi,*T,n,funct_theq(th0,w,*thi,h,soil,l,psimin),funct_theq(th1,w,*thi,h,soil,l,psimin),dfunct_theq(th0,w,*thi,h,soil,l,psimin));
			n/=2.0;
		}while(fabs(funct_theq(th0,w,*thi,h,soil,l,psimin))<fabs(funct_theq(th1,w,*thi,h,soil,l,psimin)));
		if(r==19 && c==19) printf("CR %f %f %f %f\n",th1,th0,funct_theq(th0,w,*thi,h,soil,l,psimin),dfunct_theq(th0,w,*thi,h,soil,l,psimin));
	}while(cont<maxiter && fabs(th1-th0)>tolerance);
	
	//solve equation AT+BTf(kT)+Cf(kT)+D=0 - Bisection
	if(th1!=th1 || cont==maxiter){		
		th1=soil[jres][l];
		th0=Fmax(soil[jsat][l], w/soil[jdz][l]);

		if((funct_theq(th0,w,*thi,h,soil,l,psimin)>0 && funct_theq(th1,w,*thi,h,soil,l,psimin)>0)||(funct_theq(th0,w,*thi,h,soil,l,psimin)<0 && funct_theq(th1,w,*thi,h,soil,l,psimin)<0)){
			if(r==19 && c==19) printf("\nERROR 1. IN SOIL INTERNAL ENERGY r:%ld c:%ld l:%ld\n",r,c,l);
			if(r==19 && c==19) printf("th:%30.29f theq:%f thi:%30.29f T:%f th1:%f T1:%f\n",*thw,teta_psi(Psif(*T), 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0),*thi,*T,th1,(g*(Tfreezing+tk))/(1000.0*Lf)*psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0));
			if(r==19 && c==19) printf("%f %f\n",funct_theq(th0,w,*thi,h,soil,l,psimin),funct_theq(th1,w,*thi,h,soil,l,psimin));
			if(r==19 && c==19) f=fopen(error_file_name,"a");
			if(r==19 && c==19) fprintf(f,"\nERROR 1. IN SOIL INTERNAL ENERGY r:%ld c:%ld l:%ld\n",r,c,l);
			if(r==19 && c==19) fprintf(f,"th:%30.29f thi:%30.29f T:%f th1:%f T1:%f\n",*thw,*thi,*T,th1,(g*(Tfreezing+tk))/(1000.0*Lf)*psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0));
			if(r==19 && c==19) fprintf(f,"%f %f\n",funct_theq(*thw,w,*thi,h,soil,l,psimin),funct_theq(th1,w,*thi,h,soil,l,psimin));
			if(r==19 && c==19) fclose(f);
			if(r==19 && c==19) stop_execution();
			
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
		if(r==19 && c==19) printf("\nSOL1 l:%ld h:%f thw:%f thi:%f T:%f checkE:%f\n",l,h,*thw,*thi,*T,internal_energy_soil(*thw, *thi, *T, soil[jdz][l], soil[jct][l], soil[jsat][l]));
	}else{
		*thw=th1;
		*thi=w/soil[jdz][l]-th1;
		*T=psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0)*(g*(Tfreezing+tk))/(1000.0*Lf);
		if(r==19 && c==19) printf("\nSOL2 l:%ld h:%f thw:%f thi:%f T:%f funct:%f\n",l,h,*thw,*thi,*T,funct_theq(th1,w,*thi,h,soil,l,psimin));		
		if(r==19 && c==19) printf("\n....%f %f %f %f\n",th1,soil[jres][l],psi_teta(th1, 0.0, soil[jsat][l], soil[jres][l], soil[ja][l], soil[jns][l], 1-1/soil[jns][l], psimin, 1.0),psimin);
	}
	
}*/
