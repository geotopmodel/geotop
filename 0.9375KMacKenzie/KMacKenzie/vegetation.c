
/* STATEMENT:
Vegetation module in GEOtop
This piece of code is copyrighted by the author
All rights reserved. */

//Author: Stefano Endrizzi
//Date: May 2009
//Contents: Vegetation in energy balance
#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "vegetation.h"
#include "meteo.09375.h"
#include "turbulence.h"
#include "radiation.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern long Nl, Nr, Nc;
extern double NoV;



/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void Tcanopy(long r, long c, double Tv0, double Tg, double Qg, double dQgdT, double Tg0, double Qg0, double Ta, double Qa, double zmu, double zmT,
	double z0, double z0s, double d0, double z0r, double hveg, double v, double LR, double P, double SW, double SWv, double LW, double e, double LAI,
	double *land, double Wcrn0, double Wcrnmax, double Wcsn0, double Wcsnmax, double *dWcrn, double *dWcsn, double *LWv, double *LWg, double *Hv,
	double *Hg, double *dHgdT, double *LEv, double *Eg, double *dEgdT, double *Ts, double *Qs, double *froot, double *theta, double *ftl,
	DOUBLEVECTOR *rep, PAR *par, long n, double sat, double *rh, double *rv, double *rc, double *rb, double *rh_ic, double *rv_ic, double *u_top,
	double *Etrans, double *Tv, double *Qv, double *decay, double *Locc){

	double C, C0;
	double A=1.0; /* was 0.5 see email stefano of 3/11/09*/
	double T00=Tv0, T10, T11=Tv0, T11p=Tv0, DT, Wcrn=Wcrn0, Wcsn=Wcsn0;
	double err0, err1=1.E+99, nw;
	double Lobukhov, Lobukhov0, h0=0.0, h1, dhdT, dQ;
	double subl_can, melt_can, fwliq, fwice;
	long cont, cont2, chgsgn=0;
	short a=0;

	//C0=0.02*LAI*c_liq + c_ice*Wcsn + c_liq*Wcrn;
	C0=land[jcd]*LAI*c_can + c_ice*Wcsn + c_liq*Wcrn;
	C=C0;
	fwliq=pow(Wcrn/Wcrnmax,2./3.);
	fwice=pow(Wcsn/Wcsnmax,2./3.);

	//calculates values at the instant t0 -> h0
	if(A<1){
		*Ts=0.5*Tg0+0.5*Ta;
		*Qs=0.5*Qg0+0.5*Qa;
		canopy_fluxes(r, c, T00, Tg0, Ta, Qg0, Qa, zmu, zmT, z0, z0s, d0, z0r, hveg, v, LR, P, SW, LW, e, LAI, land, Wcrn0, Wcrnmax, Wcsn0, Wcsnmax,
			&subl_can, Etrans, LWv, LWg, Hv, LEv, &h0, &dhdT, Ts, Qs, rh_ic, rv_ic, froot, theta, ftl, chgsgn, &Lobukhov, rep, par, n, sat, rh, rv,
			rc, rb, u_top, decay, Locc);
	}

	//calculates values at the instant t1 -> h1, through iterations
	cont=0;
	chgsgn=0;
	*Ts=0.5*Tg+0.5*Ta;
	*Qs=0.5*Qg+0.5*Qa;
	canopy_fluxes(r, c, T11, Tg, Ta, Qg, Qa, zmu, zmT, z0, z0s, d0, z0r, hveg, v, LR, P, SW, LW, e, LAI, land, Wcrn0, Wcrnmax, Wcsn0, Wcsnmax,
		&subl_can, Etrans, LWv, LWg, Hv, LEv, &h1, &dhdT, Ts, Qs, rh_ic, rv_ic, froot, theta, ftl, chgsgn, &Lobukhov, rep, par, n, sat, rh, rv,
		rc, rb, u_top, decay, Locc);

	melt_can = 0.0;

	do{

		T10=T11;

		//Generalized Newton-Raphson
		cont2=0;
		nw=1.0;
		err0=fabs( C*(T11-T00)/par->Dt - SWv - A*h1 - (1.-A)*h0 );

		//eq. C*(T11-T00)/Dt = cost + h(T10) + dh(T10)/dT * (T11-T10)
		//eq. C*(T11-T10)/Dt + C*(T10-T00)/Dt = cost + h(T10) + dh(T10)/dT * (T11-T10)
		//eq. T11-T10 = ( cost + h(T10) - C*(T10-T00)/Dt ) / ( C/Dt - dh/dT )

		DT=( -C*(T10-T00)/par->Dt + SWv + A*h1 + (1.0-A)*h0 ) / ( C/par->Dt - A*dhdT );

		if(DT!=DT) printf("ERROR NwRph Tcanopy T00:%f T10:%f SWv:%f h0:%f h1:%f dhdT:%f C:%f Wcsn:%f Wcrn:%f %ld %ld\n",T00,T10,SWv,h0,h1,dhdT,C,Wcsn,Wcrn,r,c);

		//printf("cont:%ld T00:%f T1:%f SWv:%f Hv:%f LEv:%f LWv:%f h0:%f h1:%f A:%f C:%f dhdt:%f DT:%f Lo:%e\n",cont,T00,T10,SWv,*Hv,*LEv,*LWv,h0,h1,A,C,dhdT,DT,Lobukhov);

		Lobukhov0=Lobukhov;

		do{

			T11 = T10 + nw*DT;

			if(subl_can<0 && T11<0){ //condensation as frost
				Wcsn=Wcsn0-subl_can*par->Dt;
				Wcrn=Wcrn0;
			}else if(subl_can<0 && T11>=0){ //condensation as dew
				Wcrn=Wcrn0-subl_can*par->Dt;
				Wcsn=Wcsn0;
			}else{	//partly evaporation, partly sublimation
				if(fwliq+fwice>0){
					Wcsn=Wcsn0-(fwice/(fwliq+fwice))*subl_can*par->Dt;
					Wcrn=Wcrn0-(fwliq/(fwliq+fwice))*subl_can*par->Dt;
				}else{
					Wcsn=Wcsn0;
					Wcrn=Wcrn0;
				}
			}

			//printf("cont2:%ld nw:%f T11:%f Wcsn:%f Wcrn:%f subl:%e\n",cont2,nw,T11,Wcsn,Wcrn,subl_can);

			if(Wcrn>Wcrnmax) Wcrn=Wcrnmax;
			if(Wcrn<0) Wcrn=0.0;
			if(Wcsn>Wcsnmax) Wcsn=Wcsnmax;
			if(Wcsn<0) Wcsn=0.0;

			if(T11>0 && Wcsn>0){	//melting
				melt_can=Fmin(Wcsn, c_ice*Wcsn*(T11-0.0)/Lf);
				T11p=T11 - Lf*melt_can/C;
				Wcsn-=melt_can;
				Wcrn+=melt_can;
			}else if(T11<0 && Wcrn>0){  //freezing
				melt_can=-Fmin(Wcrn, c_liq*Wcrn*(0.0-T11)/Lf);
				T11p=T11 - Lf*melt_can/C;
				Wcsn-=melt_can;
				Wcrn+=melt_can;
			}else{
				melt_can=0.0;
				T11p=T11;
			}

			C=land[jcd]*LAI*c_can + c_ice*Wcsn + c_liq*Wcrn;
			C=(C+C0)/2.;

			//printf("melt:%f Wcsn:%f Wcrn:%f C:%f T11p:%f\n",melt_can,Wcsn,Wcrn,C,T11p);

			canopy_fluxes(r, c, T11p, Tg, Ta, Qg, Qa, zmu, zmT, z0, z0s, d0, z0r, hveg, v, LR, P, SW, LW, e, LAI, land, Wcrn, Wcrnmax, Wcsn, Wcsnmax,
				&subl_can, Etrans, LWv, LWg, Hv, LEv, &h1, &dhdT, Ts, Qs, rh_ic, rv_ic, froot, theta, ftl, chgsgn, &Lobukhov, rep, par, n, sat, rh,
				rv, rc, rb, u_top, decay, Locc);

			err1=fabs(C*(T11-T00)/par->Dt - SWv - A*h1 - (1.0-A)*h0 );

			//if(r==59 && c==130) printf("erro:%e err1:%e Hv:%f LWv:%f LEv:%f\n",err0,err1,*Hv,*LWv,*LEv);

			nw/=3.0;
			cont2++;

		}while(err1>err0 && cont2<5);

		if(Lobukhov0*Lobukhov<0) chgsgn++;

		//if(r==59 && c==130) printf("chgsgn:%ld\n",chgsgn);

		cont++;

		if(fabs(T11-T10)<0.01 && err1<0.1) a=1;

	}while(a==0 && cont<10);

	/*if(fabs(T11-T10)>0.5){
		printf("Tcanopy not converging %f %f %ld %ld \n",T10,T11,r,c);
		if(r==59 && c==130) stop_execution();
	}*/

	*Tv=T11p;
	*dWcrn=Wcrn-Wcrn0;
	*dWcsn=Wcsn-Wcsn0;

	sat_spec_humidity(Qv, &dQ, 1.0, *Tv, P);

	turbulent_fluxes(*rh_ic, *rv_ic, P, *Ts , Tg, *Qs, Qg, dQgdT, Hg, dHgdT, Eg, dEgdT);

	if(*Tv!=(*Tv)) printf("Tv no value %ld %ld\n",r,c);

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void canopy_fluxes(long r, long c, double Tv, double Tg, double Ta, double Qg, double Qa, double zmu, double zmT, double z0, double z0s,
	double d0, double z0r, double hveg, double v, double LR, double P, double SW, double LW, double e, double LAI, double *land, double Wcrn,
	double Wcrnmax, double Wcsn, double Wcsnmax, double *Esubl, double *Etrans, double *LWv, double *LWg, double *H, double *LE, double *h,
	double *dhdT, double *Ts, double *Qs, double *rh_ic, double *rv_ic, double *froot, double *theta, double *ftl, long chgsgn, double *Lobukhov,
	DOUBLEVECTOR *rep, PAR *par, long n, double sat, double *rh, double *rv, double *rc, double *rb, double *u_top, double *decay, double *Locc){

	double rm, ft=0.0, fw, fwliq, fwice;
	double Qv, dQvdT, Hg, Lt, Lv, R, dLWvdT, dHdT, dEdT, dEsubldT, E;
	double Loc=1.E50, Loc0, Ts0;
	long cont, cont2, maxiter=5;

	//CANOPY FRACTION SET AT THE MAX OF SNOW AND LIQUID WATER FRACTION ON CANOPY
	fwliq=pow(Wcrn/Wcrnmax,2./3.);
	fwice=pow(Wcsn/Wcsnmax,2./3.);
	fw=Fmax(fwliq, fwice);
	if(fw<0) fw=0.0;
	if(fw>1) fw=1.0;

	//LONGWAVE
	longwave_vegetation(LW, e, Tg, Tv, LAI, LWv, LWg, &dLWvdT);

	//FIND SPECIFIC HUMIDITY IN THE VEGETATION LAYER
	sat_spec_humidity(&Qv, &dQvdT, 1.0, Tv, P);

	//UNDERCANOPY TURBULENT FLUXES

	//iteration for Ts
	cont2=0;
	do{

		Ts0=*Ts;

		if(chgsgn>10){	//neglects stability corrections, if the number of iterations in Tcanopy is larger than a threshold
			aero_resistance2(zmu, zmT, z0, d0, z0r, hveg, v, Ta, Ts0, Qa, *Qs, P, LR, LAI, rep, &rm, rh, rv, u_top, Lobukhov, par->state_turb, 4);
		}else{	//considers stability corrections
			aero_resistance2(zmu, zmT, z0, d0, z0r, hveg, v, Ta, Ts0, Qa, *Qs, P, LR, LAI, rep, &rm, rh, rv, u_top, Lobukhov, par->state_turb, par->monin_obukhov);
		}

		cont2++;

		//iteration for Loc (within canopy Obukhov length)
		cont=0;
		do{

			Loc0=Loc;

			if(cont==maxiter) Loc=-1.E50;

			veg_transmittance(r, c, rm, v, Ts0, Tg, z0s, LAI, z0, d0, hveg, *u_top, *Lobukhov, rb, rh_ic, rh, decay, Loc);
			*rv=*rh;

			if(Qg>(*Qs) && n==0){
				*rv_ic = (*rh_ic) + exp(8.206-4.255*sat);
			}else{
				*rv_ic = *rh_ic;
			}

			*rb=*rb/LAI;
			if(Qv<(*Qs)){	//condensation
				R=1.0;
			}else{
				canopy_evapotranspiration(*rb, Tv, Qa, P, SW, theta, land, froot, &ft, ftl);
				R=fw+(1.0-fw)*ft;
			}
			*rc=*rb/R;

			*Ts=(Ta/(*rh)+Tg/(*rh_ic)+Tv/(*rb))/(1./(*rh)+1./(*rh_ic)+1./(*rb));
			*Qs=(Qa/(*rv)+Qg/(*rv_ic)+Qv/(*rc))/(1./(*rv)+1./(*rv_ic)+1./(*rc));

			Hg=air_cp((*Ts+Tg)/2.) * air_density((*Ts+Tg)/2., (*Qs+Qg)/2., P) * (Tg-(*Ts))/(*rh_ic);

			//			-u*^3
			// Loc = ------------------    Below Canopy Monin-Obukhov length (Niu&Yang)
			//       k(g/T)(Hg/(rho*C))

			Loc=-pow(v/rm,1.5)/( ka*(g/(*Ts+tk))*(Hg/(air_density(*Ts,*Qs,P)*air_cp(*Ts))) );
			if(Hg==0.0) Loc=1.E+50;

			/*if(cont>25){
				printf("Loc:%e Lo:%e Hg:%f cont:%ld\n",Loc,*Lobukhov,Hg,cont);
			}*/

			cont++;

		}while(fabs(Loc0-Loc)>0.01 && cont<=maxiter);

		/*if(cont==maxiter){
			printf("Loc not converging, set at neutrality %ld %ld\n",r,c);
		}*/

	}while(cont2<5 && fabs((*Ts)-Ts0)>0.01);

	/*if(fabs((*Ts)-Ts0)>0.01){
		printf("Ts not converging %f %f %ld %ld\n",*Ts,Ts0,r,c);
	}*/

	//CANOPY FLUXES
	turbulent_fluxes(*rb, *rc, P, *Ts, Tv, *Qs, Qv, dQvdT, H, &dHdT, &E, &dEdT);

	//printf("H:%f E:%e v:%f cont:%ld\n",*H,E,v,cont);

	//Et from transpiration, E-Et condensation or evaporation/sublimation from water on the canopy
	Lt=Levap(Tv);
	if(E>0){	//evaporation or sublimation
		*Esubl=E*fw/R;
		dEsubldT=dEdT*fw/R;
		if(fwliq+fwice>0){
			Lv=Lt + Lf*fwice/(fwliq+fwice);	//linear interpolation to decide if sublimation or condensation occurs
		}else{
			Lv=Lt;
		}
	}else{	//condensation
		*Esubl=E;
		dEsubldT=dEdT;
		if(Tv>=0){
			Lv=Lt;
		}else{
			Lv=Lt + Lf;
		}
	}
	*Etrans=E-(*Esubl);
	*LE=Lt*(*Etrans) + Lv*(*Esubl);
	*h=(*LWv) - (*H) - (*LE);
	*dhdT=dLWvdT - dHdT - Lt*(dEdT-dEsubldT) - Lv*dEsubldT;

	if(*h!=(*h)) printf("No value in canopy fluxes Loc:%f v:%f rm:%e Ts:%f Tv:%f Ta:%f Tg:%f Hg:%f %ld %ld\n",Loc,v,rm,*Ts,Tv,Ta,Tg,Hg,r,c);

	*Locc=Loc;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void shortwave_vegetation(double Sd, double Sb, double x, double fwsn, double wsn, double Bsnd, double Bsnb, double Agd, double Agb, double C, double R, double T, double L,
	double *Sv, double *Sg){

	double phi1,phi2,G,K,w,wBd,wBb,xm;
	double b,c,d,f,h,s,u1,u2,u3,s1,s2,p1,p2,p3,p4,d1,d2,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,Iupb,Iupd,Idwb,Idwd,Ib,Id;

	phi1=0.5-0.633*C-0.33*C*C;
	phi2=0.877*(1.0-2.0*phi1);
	G=phi1+phi2*x;

	if(C==0){
		xm=1.0;
	}else{
		xm=(1.0-(phi1/phi2)*log((phi1+phi2)/phi1))/phi2;
	}

	K=G/x;

	w=(R+T)*(1.0-fwsn) + wsn*fwsn;

	wBd=(0.5*(R+T+pow(0.5*(1.0+C),2.0)*(R-T)))*(1.0-fwsn) + wsn*Bsnd*fwsn;
	wBb=( ((1.0+xm*K)/(xm*K))*0.5*w*( (G/(x*phi2+G))*(1.0-(x*phi1/(x*phi2+G))*log((x*phi1+x*phi2+G)/(x*phi1))) ) )*(1.0-fwsn) + wsn*Bsnb*fwsn;

	b=1.0-w+wBd;
	c=wBd;
	d=xm*K*wBb;
	f=w*xm*K*(1.0-wBb/w);
	h=pow(b*b-c*c,0.5)/xm;
	s=pow(xm*K,2.0)+c*c-b*b;
	u1=b-c/Agd;
	u2=b-c*Agd;
	u3=f+c*Agd;
	s1=exp(-h*L);
	s2=exp(-K*L);
	p1=b+xm*h;
	p2=b-xm*h;
	p3=b+xm*K;
	p4=b-xm*K;
	d1=p1*(u1-xm*h)/s1 - p2*(u1+xm*h)*s1;
	d2=(u2+xm*h)/s1 - (u2-xm*h)*s1;
	h1=-d*p4-c*f;
	h2=((d-h1*p3/s)*(u1-xm*h)/s1 - p2*(d-c-h1*(u1+xm*K)/s)*s2)/d1;
	h3=-((d-h1*p3/s)*(u1+xm*h)*s1 - p1*(d-c-h1*(u1+xm*K)/s)*s2)/d1;
	h4=-f*p3-c*d;
	h5=-(h4*(u2+xm*h)/(s*s1) + (u3-h4*(u2-xm*K)/s)*s2)/d2;
	h6=(h4*(u2-xm*h)*s1/s + (u3-h4*(u2-xm*K)/s)*s2)/d2;
	h7=c*(u1-xm*h)/(d1*s1);
	h8=-c*(u1+xm*h)*s1/d1;
	h9=(u2+xm*h)/(d2*s1);
	h10=-s1*(u2-xm*h)/d2;

	Iupb=h1/s+h2+h3;
	Iupd=h7+h8;
	Idwb=h4*exp(-K*L)/s+h5*s1+h6/s1;
	Idwd=h9*s1+h10/s1;

	Ib=(1.0-Iupb) - (1.0-Agd)*Idwb - (1.0-Agb)*exp(-K*L);
	Id=(1.0-Iupd) - (1.0-Agd)*Idwd;

	if(x>0){
		*Sv=Sb*Ib + Sd*Id;
		*Sg=Sb*(1.0-Agb)*exp(-K*L) + (Sb*Idwb+Sd*Idwd)*(1.0-Agd);
	}else{
		*Sv=Sd*Id;
		*Sg=Sd*Idwd*(1.0-Agd);
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void longwave_vegetation(double Lin, double eg, double Tg, double Tv, double L, double *Lv, double *Lg, double *dLv_dTv){

	double ev, Lvdw;

	ev=1.0-exp(-L);

	Lvdw=(1.0-ev)*Lin + ev*SB(Tv);

	*Lv=ev*(1.0+(1.0-eg)*(1.0-ev))*Lin + ev*eg*SB(Tg) - (2.0-ev*(1.0-eg))*ev*SB(Tv);
	*Lg=eg*Lvdw - eg*SB(Tg);
	*dLv_dTv=-(2.0-ev*(1.0-eg))*ev*dSB_dT(Tv);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void canopy_rain_interception(double rain_max_loading, double LAI, double Prain, double *max_storage, double *storage, double *drip){

	*max_storage = 0.1*LAI; //mm
	*storage = (*storage) + Prain*(1.0-exp(-0.5*LAI));
	*drip=Prain*exp(-0.5*LAI);

	if( (*storage) > (*max_storage) ){
		*drip = *drip + (*storage) - (*max_storage);
		*storage = (*max_storage);
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void canopy_snow_interception(double snow_max_loading, double LAI, double Psnow, double Tc, double v, double Dt, double *max_storage, double *storage, double *drip){

	double load, unload;
	double CT=1.87E5, CV=1.56E5;

	*max_storage = snow_max_loading*LAI;

	load = ( (*max_storage) - (*storage) ) * (1. - exp(-Psnow/(*max_storage)));	//Niu & Yang, 2004
	*storage = (*storage) + load;
	*drip = Psnow - load;

	if((*storage)>(*max_storage)){
		*drip = *drip - ( (*max_storage) - (*storage) );
		*storage = (*max_storage);
	}

	unload = (*storage)*(Fmax(0.0, Tc+3.0)/CT + v/CV)*Dt;
	if(unload<0.1) unload=0.0;	//prevents very low snowfalls

	*drip = *drip + unload;

	*storage = (*storage) - unload;

}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void update_roughness_veg(double hc, double snowD, double zmu, double zmt, double *z0_ris, double *d0_ris, double *hc_ris){


	*hc_ris=hc-1.E-3*snowD;
	*d0_ris=0.667*(*hc_ris);
	//*hc_ris=hc;
	//*d0_ris=0.0;
	*z0_ris=0.1*(*hc_ris);

	if(zmu<*hc_ris) t_error("zmu lower than hc");
	if(zmt<*hc_ris) t_error("zmt lower than hc");
	if(*hc_ris-(*d0_ris)<(*z0_ris)) t_error("hc-d0 lower than z0");
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void root(double d, double *D, double *root_fraction){

	long l;
	double z=0.0;

	for(l=1;l<=Nl;l++){
		z+=D[l];
		if(d>z){
			root_fraction[l]=D[l]/d;
		}else{
			if(d>(z-D[l])){
				root_fraction[l]=(d-(z-D[l]))/d;
			}else{
				root_fraction[l]=0.0;
			}
		}
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void canopy_evapotranspiration(double rbv, double Tv, double Qa, double Pa, double SWin, double *theta, double *land, double *root, double *f, double *fl){

	double fS, fe, fTemp, Rsmin, ea, ev, de;
	long l;

	//CANOPY TRANSPIRATION (parameters from Best (1998))
	//solar radiation [Best, (1998); Dolman et al., 1991]
	fS=SWin/(SWin+250.0)*1.25;

	//pressure deficit [Best, (1998); Dickinson et al., 1991]*/
	ea=Qa*Pa/(0.378*Qa+0.622);
	sat_vap_pressure(&ev, &de, Tv, Pa);
	fe=1.0+(ev-ea)/40.0;

	//temperature [Best, (1998); Dickinson et al., 1991]*/
	if(Tv<=0){
		fTemp=1E-12;
	}else if(Tv>=50.0){
		fTemp=1E-12;
	}else{
		fTemp=(Tv-0.0)*(50.0-Tv)/625.0;
	}

	*f=0.0;
	for(l=1;l<=Nl;l++){
		//water content [Wigmosta et al., (1994); Feddes et al.(1978)]
		if (theta[l]>=land[jtfc]){
			fl[l]=1.0;
		}else if(theta[l]>land[jtwp]){
			fl[l]=(theta[l]-land[jtwp])/(land[jtfc]-land[jtwp]);
		}else{
			fl[l]=0.0;
		}

		//stomata resistance for each layer
		if(fS*fe*fTemp*fl[l]<6.0E-11){
			fl[l]=1.0E12;
		}else{
			Rsmin=land[jrs];
			fl[l]=Rsmin/(fS*fe*fTemp*fl[l]);
		}

		//transpiration fraction for each layer
		fl[l]=root[l]*(rbv/(rbv+fl[l]));

		//transpiration for all the column (as a fraction of Epc)
		*f+=fl[l];
	}

	for(l=1;l<=Nl;l++){
		if(*f!=0) fl[l]/=(*f);
	}

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void veg_transmittance(long r, long c, double rm, double v, double Ts, double Tg, double z0soil, double LAI, double z0veg, double d0veg,
	double Hveg, double u_top, double Lo, double *rb, double *rh, double *ra, double *decay, double Loc){

	//double Cs, Cs_bare, Cs_dense = 0.004;
	double Lc = 0.4;	//characteristic dimension of vegetation [m]
	double u_star = pow(v/rm,0.5);
	double u_veg;
	double zm = d0veg + z0veg;
	//double W = 1.0-exp(-LAI);
	double phi=1.0, philog;
	//double W = 1.0;
	double n = ndc;

	//stability
	//over canopy
	if(Lo<0){
		phi=pow(1.0-16.0*Hveg/Lo,-0.5);
	}else{
		phi=1.0+5.0*Hveg/Lo;
	}
	phi=1.0;

	//under canopy
	if(Loc<0){
		philog=pow(1.0-15.0*zm/Loc,-0.25);
	}else{
		philog=1.0+4.7*zm/Loc;
	}
	n=Fmin(1.E+5, n*pow(philog,0.5));
	*decay=n;

	//update aerodynamic resistance
	//*ra=*ra + (Hveg/(n*ka*u_star*(Hveg-d0veg)/phi))*(exp(n*(1.0-zm/Hveg))-1.0);

	//Zeng
	//u_veg=u_star;

	//Huntingford
	u_veg=Fmax(0.001,u_top*exp(n*(zm/Hveg-1.0)));

	//canopy resistance
	*rb=Fmin(1.E20, 70.0*pow(Lc/u_veg,0.5));

	//ground resistance
	//Zeng
	//Cs_bare=(pow(z0soil*u_star/1.5E-5,-0.45)*ka/0.13);
	//*rh=1.0/( (W*Cs_dense + (1.0-W)*Cs_bare) * u_star);

	//Huntingford
	*rh=Fmin(1.E20,(Hveg*exp(n)/(n*ka*u_star*(Hveg-d0veg)/phi))*(exp(-n*z0soil/Hveg)-exp(-n*zm/Hveg)));
	if(n<10) *rh=*rh+pow((exp(n)-exp(n*(1.-d0veg/Hveg)))*Hveg/(n*d0veg),0.45)*(u_star*2./ka);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
