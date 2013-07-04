#include "keywords_file.h"
#include "constant.h"
#include "struct.geotop.09375.h"
#include "turbulence.h"
#include "meteo.09375.h"

extern T_INIT *UV;
extern char *WORKING_DIRECTORY;
extern STRINGBIN *files;
extern char *error_file_name;
extern long Nl, Nr, Nc;
extern double NoV;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, double T, double Qa, double Q, double P, double gmT, DOUBLEVECTOR *rep, double *rm, double *rh, 
	double *rv, short state_turb, short MO){

	//calculate resistences
	
	//double p=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	double p=1.0;
	double Tpa=(Ta+tk)*p,Tp=(T+tk)*p;	//potential temperatures

	if(state_turb==0){
		Lewis(zmu, zmt, d0, z0, z0_z0t, Tpa, Tp, v, rm, rh, rv, rep);
	}else if(state_turb==1){
		Businger(MO, zmu, zmt, d0, z0, v, 0.5*T+0.5*Ta, T-Ta, Q-Qa, z0_z0t, rm, rh, rv, rep);		
	}else if(state_turb==2){	//catabatic flows - OERLEMANS & GRISOGONO (2002)
		*rm=1.0/(v*ka*ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 )));
		if(p*(-gmT+0.0098)>0.0015){
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*p*(-gmT+0.0098)*5),0.5) );
		}else{
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*(0.0015)*5),0.5) );
		}
		*rv=*rh;
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void aero_resistance2(double zmu, double zmt, double z0, double d0, double z0_z0t, double hveg, double v, double Ta, double T, double Qa, double Q, double P,
	double gmT, double LAI, DOUBLEVECTOR *rep, double *rm, double *rh, double *rv, double *u_top, double *Lmo, short state_turb, short MO){

	//calculate resistences
	
	//double p=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	double p=1.0;
	double Tpa=(Ta+tk)*p,Tp=(T+tk)*p;	//potential temperatures
	double u_star;
//	double n = ndc;

	if(state_turb==0){
		Lewis(zmu, zmt, d0, z0, z0_z0t, Tpa, Tp, v, rm, rh, rv, rep);
		u_star=pow(v/(*rm),0.5);
		*u_top=(u_star/ka)*log((hveg-d0)/z0);
		
	}else if(state_turb==1){
		Businger(MO, zmu, zmt, d0, z0, v, 0.5*T+0.5*Ta, T-Ta, Q-Qa, z0_z0t, rm, rh, rv, rep);	
		*Lmo=rep->co[2];
		u_star=pow(v/(*rm),0.5);
		*u_top=(u_star/ka)*CZ(MO, hveg, z0, d0, *Lmo, (*Psim));

	}else if(state_turb==2){	//catabatic flows - OERLEMANS & GRISOGONO (2002)
		*rm=1.0/(v*ka*ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 )));
		if(p*(-gmT+0.0098)>0.0015){
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*p*(-gmT+0.0098)*5),0.5) );
		}else{
			*rh=1/( 0.0004*(Ta-T)*pow(g/(tk*(0.0015)*5),0.5) );
		}
		*rv=(*rh);
		*u_top=v; //fictitious value, as it is not used
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void turbulent_fluxes(double rh, double rv, double P, double Ta, double T, double Qa, double Q, double dQdT, double *H, double *dHdT, double *E, double *dEdT){

	double rho, cp, pot;

	rho=air_density(0.5*(Ta+T), 0.5*(Qa+Q), P);
	cp=air_cp(0.5*(Ta+T));

	//pot=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	pot=1.0;

	//sensible heat flux [W/m2]
	/*a maximum value of the resistance for the sensible heat flux was introduced according to Jordan et el., 1999
	(HEAT BUDGET OF SNOW-COVERED SEA ICE AT NORTH POLE 4) as Cwindless=0.5 W m^2 K^-1, rh=rho*cp/C=1300/0.5=2600*/

	//*H=cp*rho*pot*(T-Ta)/Fmin(rh,2.6E3);
	//*dHdT=cp*rho*pot/Fmin(rh,2.6E3);
	*H=cp*rho*pot*(T-Ta)/rh;
	*dHdT=cp*rho*pot/rh;
	//evaporation [kg/(s*m2)]
	*E=rho*(Q-Qa)/rv;
	*dEdT=rho*dQdT/rv;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//****PSIm
double Psim(double z)
{
double x,psi;
x=pow(1.0-15.0*z,0.25);
psi=2.0*log((1.0+x)/2.0)+log((1.0+x*x)/2.0)-2.0*atan(x)+0.5*Pi;
return(psi);
}

//****Psih
double Psih(double z)
{
double x,psi;
x=pow(1.0-15.0*z,0.25);
psi=2.0*log((1.0+x*x)/2.0);
return(psi);
}

//****Zero
double Zero(double z)
{
return(0.0);
}

//****PsiHolstag&deBruin
double PsiStab(double z)
{
double psi;
//psi=5*z;
psi=10.71 + 0.7*z + 0.75*(z-14.28)*exp(-0.35*z);	//Holstag&De Bruin
/*if(z<=1){	//Brutsaert
	psi=5*z;
}else{
	psi=5*(1+log(z));
}*/
return(psi);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//Loius' scheme (Kot & Song, 1998)
void Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double z0t, f, Rib;
	double Chn, bh, c1h, c2h, c3h, c4h, Chx, ch, Fh;
	double Cmn, bm, c1m, c2m, c3m, c4m, Cmx, cm, Fm;

	//check
	//if(zmt-zmu>0.5 || zmt-zmu<-0.5) t_error("If you use Louis' scheme, wind and temperature should be measured approximately at the same elevation on the ground");

	//roughness
	if(z0_z0t==0.0){	//rigid surface
		z0t=z0/8.5;
	}else{				//bending surface
		z0t=z0/z0_z0t;
	}

	/* CH neutrale [m/s]*/
	Cmn=ka*ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 ));
	Chn=ka*ka/(log((zmu-d0)/z0 )*log((zmt-d0)/z0t));

	/* calcola la funzione di stabilita' secondo la trattazione semplificata di Garrat, 1992 */
	f=(zmu/zmt-z0/zmt)/pow(zmt-z0t,0.5);

	/*calculation of the Richardson'number of Bulk */
	Rib=9.81*zmt*(Ta-Ts)/(0.5*(Ta+Ts)*v*v)*f*f;

	bm=8.0;
	bh=23.0;
	if (Rib < 0.0) {
		c1m=-0.9848;
		c2m=2.5398;
		c3m=-0.2325;
		c4m=14.1727;
		if (z0/z0t<100.0) {
			c1h= -1.1790;
			c2h= -1.9256;
			c3h= 0.1007;
			c4h= 16.6796;
		}else{
			c1h= -1.0487;
			c2h= -1.0689;
			c3h= 0.0952;
			c4h= 11.7828;
		}
	}else{
		if (z0/z0t<100.0) {
			c1m=-0.4738;
			c2m=-0.3268;
			c3m=0.0204;
			c4m=10.0715;
			c1h= -0.5128;
			c2h= -0.9448;
			c3h= 0.0643;
			c4h= 10.8925;
		}else{
			c1m=-0.4613;
			c2m=-0.2402;
			c3m=0.0146;
			c4m=8.9172;
			c1h= -0.3169;
			c2h= -0.3803;
			c3h= 0.0205;
			c4h= 7.5213;
		}
	}

	Cmx=c1m*log(zmu/z0)+c2m*log(z0/z0t)+c3m*log(zmt/z0t)+c4m;
	Chx=c1h*log(zmu/z0)+c2h*log(z0/z0t)+c3h*log(zmt/z0t)+c4h;

	cm=Cmx*Cmn*bm*f*pow(pow(zmu/z0 ,1.0/3.0)-1.0,1.5);
	ch=Chx*Chn*bh*f*pow(pow(zmt/z0t,1.0/3.0)-1.0,1.5);

	if (Rib < 0.0) {
		Fm=1.0-bm*Rib/(1.0+cm*pow((-Rib),0.5));
		Fh=1.0-bh*Rib/(1.0+ch*pow((-Rib),0.5));
	}else{
		Fm=1.0/(pow(1.0+Cmx*Rib,2.0));
		Fh=1.0/(pow(1.0+Chx*Rib,2.0));
	}

	/* calcola la resistenza aereodinamica */
	*rm=1.0/(Cmn*Fm*v);
	*rh=1.0/(Chn*Fh*v);
	*rv=*rh;

	//check variables
	w->co[1]=Rib;
	w->co[2]=Fh;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double cz(double zmeas, double z0, double d0, double L, double (* unstab)(double z), double (* stab)(double z)){

	double c,zeta;

	zeta=(zmeas-d0)/L;

	if(zeta<0){
		c=log((zmeas-d0)/z0) - (*unstab)(zeta) + (*unstab)(z0/L);
	}else{
		c=log((zmeas-d0)/z0) + (*stab)(zeta) - (*stab)(z0/L);
	}

	return(c);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(double z)){

	double c=0; /* modified by Emanuele Cordano on 24/9/9 */

	if(state==1){			//both instability and stability considered
		c=cz(zmeas,z0,d0,L,(Psi),(*PsiStab));
	}else if(state==2){		//instability considered & stability not considered
		c=cz(zmeas,z0,d0,L,(Psi),(*Zero));
	}else if(state==3){		//instability not considered & stability considered
		c=cz(zmeas,z0,d0,L,(*Zero),(*PsiStab));
	}else if(state==4){		//both instability and stability not considered
		c=cz(zmeas,z0,d0,L,(*Zero),(*Zero));
	}else{
		t_error("Value not admitted in CV");
	}

	return(c);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, double M, double N, double R, double *var, double *c, double *z0v,
	double (*Psi)(double z), double (*roughness)(double x, double y, double z) ){

	*z0v=z0*(*roughness)(M, N, R);
	//if(*z0v<1.0E-5) *z0v=1.0E-5;
	*c=CZ(a,zmeas,*z0v,d0,L,(Psi));
	*var=delta*ka/(*c);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double roughT(double M, double N, double R){

	double b0,b1,b2,fr;

	if(M<=0.135){
		b0=1.250;
		b1=0.0;
		b2=0.0;
	}else if(M<2.5){
		b0=0.149;
		b1=-0.550;
		b2=0.0;
	}else{
		b0=0.317;
		b1=-0.565;
		b2=-0.183;
	}
	fr=R+N*exp(b0+b1*log(M)+b2*pow(log(M),2.0));

	return(fr);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double roughQ(double M, double N, double R){

	double b0,b1,b2,fr;

	if(M<=0.135){
		b0=1.610;
		b1=0.0;
		b2=0.0;
	}else if(M<2.5){
		b0=0.351;
		b1=-0.628;
		b2=0.0;
	}else{
		b0=0.396;
		b1=-0.512;
		b2=-0.180;
	}
	fr=R+N*exp(b0+b1*log(M)+b2*pow(log(M),2.0));
	//fr=R+N*exp(-ka*(7.3*pow(M,0.25)*pow(0.595,0.5)-5));
	return(fr);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double L, L0, cm, u_star, ch, T_star, cv, Q_star;
	double z0v, z0t, z0q;
	long cont;
	double tol;	
	FILE *f;
	
	//first guess of Obukhov length
	if(DT<0){
		L=1.E5;
	}else{
		L=-1.E5;
	}

	cont=0;
	tol=1.0E99;

	do{
	
		if(cont>0) tol=10*T_star+100*u_star+1000*Q_star;	

		L0=L;

		//Conductances
		Star(a, zmu, z0, d0, L, 0.0, v, 1.0, 0.0, 1.0, &u_star, &cm, &z0v, (*Psim), (*roughT));	//momentum
		if(z0_z0t==0.0){	//rigid surface
			Star(a, zmt, z0, d0, L, u_star, DT, u_star*z0/1.4E-5, 1.0, 0.0, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, u_star*z0/1.4E-5, 1.0, 0.0, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}else{	//bending surface
			Star(a, zmt, z0, d0, L, u_star, DT, 1.0, 0.0, 1.0/z0_z0t, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, 1.0, 0.0, 1.0/z0_z0t, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}

		//Obukhov length
		L=-u_star*u_star*(T+tk)/(ka*g*(T_star-0.61*Q_star*(T+tk)));
		
		//printf("v:%f u_star:%f T_star:%f DT:%f L:%e\n",v,u_star,T_star,DT,L);

		cont++;

	}while(fabs(100*T_star+100*u_star+1000*Q_star-tol)>0.01 && cont<=10);

	if(d0>zmu || d0>zmt){
		f=fopen(error_file_name,"a");
		fprintf(f,"ERROR: Displacement height greater than measurement elevations");
		fclose(f);
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f=fopen(error_file_name,"a");
		fprintf(f,"ERROR: Elevation of sensors lower than roughness length: zmu=%f zmt=%f z0=%f z0t=%f\n",zmu,zmt,z0,z0t);
		fclose(f);
	}

	*rm=cm*cm/(ka*ka*v);
	*rh=ch*cm/(ka*ka*v);
	*rv=cv*cm/(ka*ka*v);

	w->co[1]=(double)cont;
	w->co[2]=L;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Businger2(long r, long c, short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, double *rm, double *rh, double *rv, DOUBLEVECTOR *w){

	double L, L0, cm, u_star, ch, T_star, cv, Q_star;
	double z0v, z0t, z0q;
	double tol;
	long cont;
	FILE *f;

	//first guess of Obukhov length
	if(DT<0){
		L=1.E5;
	}else{
		L=-1.E5;
	}

	cont=0;
	tol=1.0E99;

	do{

		if(cont>0) tol=10*T_star+100*u_star+1000*Q_star;

		L0=L;

		//Conductances
		Star(a, zmu, z0, d0, L, 0.0, v, 1.0, 0.0, 1.0, &u_star, &cm, &z0v, (*Psim), (*roughT));	//momentum
		if(z0_z0t==0.0){	//rigid surface
			Star(a, zmt, z0, d0, L, u_star, DT, u_star*z0/1.4E-5, 1.0, 0.0, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, u_star*z0/1.4E-5, 1.0, 0.0, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}else{	//bending surface
			Star(a, zmt, z0, d0, L, u_star, DT, 1.0, 0.0, 1.0/z0_z0t, &T_star, &ch, &z0t, (*Psih), (*roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, 1.0, 0.0, 1.0/z0_z0t, &Q_star, &cv, &z0q, (*Psih), (*roughQ)); //water vapour flux
		}

		//Obukhov length
		L=-u_star*u_star*(T+tk)/(ka*g*(T_star+0.61*Q_star*(T+tk)));
		if(a==4) L=1.E5;
		//if(L*L0<0) L=L0;


		cont+=1;

	}while(fabs(10*T_star+100*u_star+1000*Q_star-tol)>0.01 && cont<=100);

	if(d0>zmu || d0>zmt){
		f=fopen(error_file_name,"a");
		fprintf(f,"ERROR: Displacement height greater than measurement elevations");
		fclose(f);
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f=fopen(error_file_name,"a");
		fprintf(f,"ERROR: Elevation of sensors lower than roughness length: zmu=%f zmt=%f z0=%f z0t=%f\n",zmu,zmt,z0,z0t);
		fclose(f);
	}

	*rm=cm*cm/(ka*ka*v);
	*rh=ch*cm/(ka*ka*v);
	*rv=cv*cm/(ka*ka*v);

	w->co[1]=(double)cont;
	w->co[2]=L;

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Levap(double T){
	double Lv;
	if(T>0.0){
		Lv=2501000.0+(2406000.0-2501000.0)/40.0*T;
	}else{
		Lv=2501000.0;
	}
	return(Lv);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double latent(double Ts, double Le){
	double L;
	if(Ts<0){
		L=Le+Lf;	//sublimazione
	}else{
		L=Le;		//evaporazione
	}
	return(L);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
