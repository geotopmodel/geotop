
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
#include "turbulence.h"
#include "geotop_common.h"

using namespace std;

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Turbulence::aero_resistance(double zmu, double zmt, double z0, double d0, double z0_z0t, double v, double Ta, double T, double Qa, 
						   double Q, double P, double gmT, double *Lobukhov, double *rm, double *rh, double *rv, short state_turb, 
						   short MO, long maxiter) 
{
	FILE *f;
	
	//calculates resistences
	
	//double p=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	double p=1.0;
	//double Tpa=(Ta+tk)*p,Tp=(T+tk)*p;	//potential temperatures

	if(state_turb==0){
		//Lewis(zmu, zmt, d0, z0, z0_z0t, Tpa, Tp, v, rm, rh, rv, rep);
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: state_turb == 0 not possible, check option file\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
		
	}else if(state_turb==1){
		Businger(MO, zmu, zmt, d0, z0, v, 0.5*T+0.5*Ta, T-Ta, Q-Qa, z0_z0t, rm, rh, rv, Lobukhov, maxiter);	
		
	}else if(state_turb==2){	//catabatic flows - OERLEMANS & GRISOGONO (2002)
		*rm=1.0/(v*GTConst::ka*GTConst::ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 )));
		if(p*(-gmT+0.0098)>0.0015){
			*rh=1/( 0.0004*(Ta-T)*pow(GTConst::GRAVITY/(GTConst::tk*p*(gmT+0.0098)*5),0.5) );
		}else{
			*rh=1/( 0.0004*(Ta-T)*pow(GTConst::GRAVITY/(GTConst::tk*(0.0015)*5),0.5) );
		}
		*rv=*rh;
	}
	
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Turbulence::turbulent_fluxes(const double& rh, const double& rv, const double& P, const double& Ta, const double& T, const double& Qa,
                                  const double& Q, const double& dQdT, double& H, double& dHdT, double& E, double& dEdT)
{
	double rho, cp, pot;
	
	rho = air_density(0.5*(Ta+T), Qa, P);
	cp  = air_cp(0.5*(Ta+T));
	
	//pot=pow((1000.0/P),(0.286*(1-0.23*Qa)));
	pot = 1.0;
	
	//sensible heat flux [W/m2] 
	/*a maximum value of the resistance for the sensible heat flux was introduced according to Jordan et el., 1999
	(HEAT BUDGET OF SNOW-COVERED SEA ICE AT NORTH POLE 4) as Cwindless=0.5 W m^2 K^-1, rh=rho*cp/C=1300/0.5=2600*/
	H    = cp * rho * pot * (T-Ta) / Fmin(rh, 2.6E3);	
	dHdT = cp * rho * pot / Fmin(rh, 2.6E3);	
	//*H=cp*rho*pot*(T-Ta)/rh;	
	//*dHdT=cp*rho*pot/rh;		
	
	//evaporation [kg/(s*m2)]
	E    = rho * (Q-Qa) /rv;	
	dEdT = rho * dQdT /rv;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//****PSIm
double Turbulence::Psim(const double& z)
{
	double x,psi;
	x=pow(1.0-15.0*z,0.25);
	psi=2.0*log((1.0+x)/2.0)+log((1.0+x*x)/2.0)-2.0*atan(x)+0.5*GTConst::Pi;

	return psi;
}

//****Psih
double Turbulence::Psih(const double& z)
{
	double x,psi;
	x=pow(1.0-15.0*z,0.25);
	psi=2.0*log((1.0+x*x)/2.0);
	return psi;
}

//****Zero
double Turbulence::Zero(const double&)
{
	return 0.0;
}

//****PsiHolstag&deBruin
double Turbulence::PsiStab(const double& z)
{
	double psi = 10.71 + 0.7*z + 0.75*(z-14.28)*exp(-0.35*z);	//Holstag&De Bruin

	return psi;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

//Loius' scheme (Kot & Song, 1998)
//void Turbulence::Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh,
//				   double *rv, DOUBLEVECTOR *w)
void Turbulence::Lewis(double zmu, double zmt, double d0, double z0, double z0_z0t, double Ta, double Ts, double v, double *rm, double *rh, 
				   double *rv, GeoVector<double>& w)
{
	double z0t, f, Rib;
	double Chn, bh, c1h, c2h, c3h, c4h, Chx, ch, Fh;
	double Cmn, bm, c1m, c2m, c3m, c4m, Cmx, cm, Fm;

	//roughness
	if(z0_z0t==0.0){	//rigid surface
		z0t=z0/8.5;
	}else{				//bending surface
		z0t=z0/z0_z0t;
	}

	/* CH neutrale [m/s]*/
	Cmn=GTConst::ka*GTConst::ka/(log((zmu-d0)/z0 )*log((zmu-d0)/z0 ));
	Chn=GTConst::ka*GTConst::ka/(log((zmu-d0)/z0 )*log((zmt-d0)/z0t));

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
	w[1]=Rib;
	w[2]=Fh;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Turbulence::cz(double zmeas, double z0, double d0, double L, double (* unstab)(const double& z), double (* stab)(const double& z))
{
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

double Turbulence::CZ(short state, double zmeas, double z0, double d0, double L, double (*Psi)(const double& z))
{
	double c;
	FILE *f;

	if(state==1){			//both instability and stability considered
		c=cz(zmeas,z0,d0,L,(Psi),(&Turbulence::PsiStab));
	}else if(state==2){		//instability considered & stability not considered
		c=cz(zmeas,z0,d0,L,(Psi),(&Turbulence::Zero));
	}else if(state==3){		//instability not considered & stability considered
		c=cz(zmeas,z0,d0,L,(*Zero),(&Turbulence::PsiStab));
	}else if(state==4){		//both instability and stability not considered
		c=cz(zmeas,z0,d0,L,(*Zero),(&Turbulence::Zero));
	}else{
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Value of state turbulence not admitted\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
		c = 0.0;
	}

	return(c);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Turbulence::Star(short a, double zmeas, double z0, double d0, double L, double u, double delta, double M, double N, double R, 
		  double *var, double *c, double *z0v, double (*Psi)(const double& z), double (*roughness)(const double& x, const double& y, const double& z))
{	
	*z0v=z0*(*roughness)(M, N, R);
	//if(*z0v<1.0E-5) *z0v=1.0E-5;
	*c=CZ(a,zmeas,*z0v,d0,L,(Psi));
	*var=delta*GTConst::ka/(*c);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Turbulence::roughT(const double& M, const double& N, const double& R)
{
	double b0,b1,b2,fr;

	if (M<=0.135) {
		b0=1.250;
		b1=0.0;
		b2=0.0;
	} else if (M<2.5) {
		b0=0.149;
		b1=-0.550;
		b2=0.0;
	} else {
		b0=0.317;
		b1=-0.565;
		b2=-0.183;
	}

	fr = R + N * exp(b0+b1*log(M)+b2*pow(log(M),2.0));
	
	return(fr);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Turbulence::roughQ(const double& M, const double& N, const double& R)
{
	double b0,b1,b2,fr;
	
	if (M<=0.135) {
		b0=1.610;
		b1=0.0;
		b2=0.0;
	} else if (M<2.5) {
		b0=0.351;
		b1=-0.628;
		b2=0.0;
	} else {
		b0=0.396;
		b1=-0.512;
		b2=-0.180;
	}

	fr = R + N * exp(b0+b1*log(M)+b2*pow(log(M),2.0));
	//fr=R+N*exp(-ka*(7.3*pow(M,0.25)*pow(0.595,0.5)-5));
	return(fr);
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Turbulence::Businger(short a, double zmu, double zmt, double d0, double z0, double v, double T, double DT, double DQ, double z0_z0t, 
			  double *rm, double *rh, double *rv, double *Lobukhov, long maxiter)
{

	double L, cm, u_star, ch, T_star, cv, Q_star;
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
		
		//Conductances
		Star(a, zmu, z0, d0, L, 0.0, v, 1.0, 0.0, 1.0, &u_star, &cm, &z0v, (&Turbulence::Psim), (&Turbulence::roughT));	//momentum
		if(z0_z0t==0.0){	//rigid surface
			Star(a, zmt, z0, d0, L, u_star, DT, u_star*z0/1.4E-5, 1.0, 0.0, &T_star, &ch, &z0t, (*Psih), (&Turbulence::roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, u_star*z0/1.4E-5, 1.0, 0.0, &Q_star, &cv, &z0q, (*Psih), (&Turbulence::roughQ)); //water vapour flux
		}else{	//bending surface
			Star(a, zmt, z0, d0, L, u_star, DT, 1.0, 0.0, 1.0/z0_z0t, &T_star, &ch, &z0t, (*Psih), (&Turbulence::roughT)); //heat flux
			Star(a, zmt, z0, d0, L, u_star, DQ, 1.0, 0.0, 1.0/z0_z0t, &Q_star, &cv, &z0q, (*Psih), (&Turbulence::roughQ)); //water vapour flux
		}

		//Obukhov length
		L=-u_star*u_star*(T+GTConst::tk)/(GTConst::ka*GTConst::GRAVITY*(T_star-0.61*Q_star*(T+GTConst::tk)));
				
		cont++;

	}while(fabs(100*T_star+100*u_star+1000*Q_star-tol)>0.01 && cont<=maxiter);
	
	if(d0>zmu || d0>zmt){
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Displacement height greater than measurement elevations\n");
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}
	if(zmu<=z0 || zmt<=z0t || zmt<=z0q){
		f = fopen(geotop::common::Variables::FailedRunFile.c_str(), "w");
		fprintf(f,"Error:: Elevation of sensors lower than roughness length: zmu=%f zmt=%f z0=%f z0t=%f\n",zmu,zmt,z0,z0t);
		fclose(f);
		t_error("Fatal Error! Geotop is closed. See failing report.");	
	}

	*rm=cm*cm/(GTConst::ka*GTConst::ka*v);
	*rh=ch*cm/(GTConst::ka*GTConst::ka*v);
	*rv=cv*cm/(GTConst::ka*GTConst::ka*v);
	
	*Lobukhov=L;
	
}


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Turbulence::Levap(const double& T)
{
	double Lv;

	if (T > 0.0) {
		Lv = 2501000.0+(2406000.0-2501000.0) / 40.0 * T;
	}else{
		Lv = 2501000.0;
	}

	return Lv;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

double Turbulence::latent(const double& Ts, const double& Le)
{
	double L;

	if (Ts<0) {
		L = Le + GTConst::Lf;	//sublimation
	} else {
		L = Le;        		//evaporation
	}

	return L;
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void Turbulence::find_actual_evaporation_parameters(long R, long C, double *alpha, double *beta, GeoVector<double>& evap_layer, double *theta,
										  //double **soil, double *T, double psi, double P, double rv, double Ta, double Qa, 
										  GeoTensor<double>& soil, long sy, double *T, double psi, double P, double rv, double Ta, double Qa, 
										  double Qgsat, long nsnow)
{	
//	DOUBLEVECTOR *r;
	GeoVector<double> r;
	double hs, F, A, B, Qs, D, Qsat, rho;
//	long l, n = evap_layer->nh;
	long l, n = evap_layer.size();
	
	//fprintf(stderr, "THOEMS: evap_layer.size()=%d\n" , evap_layer.size());
	
	if(nsnow>0){//snow or ice
		
		*alpha = 1.0;
		*beta = 1.0;
	//	initialize_doublevector(evap_layer, 0.);
		evap_layer.resize(evap_layer.size(),0.0);
				
	}else{
		
		rho = air_density(0.5*(Ta+T[1]), Qa, P);
		
		//from Ye and Pielke, 1993 - bare soil evaporation
		
		if(psi > 10.){		//ponding
			
			*alpha = 1.0;
			*beta = 1.0;
			
		//	evap_layer->co[1] = rho * ( Qgsat - Qa ) / rv;
			evap_layer[1] = rho * ( Qgsat - Qa ) / rv;
			
		}else if(theta[1] >= soil[sy][jsat][1]){	//saturation
			
			*alpha = 1.0;
			*beta = theta[1];
			
		//	evap_layer->co[1] = theta[1] * rho * ( Qgsat - Qa ) / rv;
			evap_layer[1] = theta[1] * rho * ( Qgsat - Qa ) / rv;

		}else{	//unsaturated
			
			A = 0.;
			B = 0.;
			
			//r = new_doublevector(n);	////molecular diffusivity resistances [s/m]
			r.resize(n+1);

			//calculates water vapor fluxes
			for(l=1;l<n;l++){
											
				D = GTConst::D00 * pow((T[l]+GTConst::tk)/GTConst::tk, 2.) * (GTConst::Pa0/P);	//molecular diffusivity water vapor [mm2/s]
				r[l]  = (1.E3/D) * soil[sy][jdz][l];  
				//r->co[l] = (1.E3/D) * soil[jdz][l];  
				
				Qsat = SpecHumidity(SatVapPressure(T[l], P), P);
				//if(l>1) r->co[l] += r->co[l-1];
				if (l > 1) r[l] += r[l-1];
			
				if( theta[l] <= soil[sy][jfc][l] ){
					hs = 0.5 * ( 1. - cos( GTConst::Pi * ( theta[l] - soil[sy][jres][l] ) / ( soil[sy][jfc][l] - soil[sy][jres][l] ) ) );
				}else{
					hs = 1.;
				}			

				//evap_layer->co[l] = rho*(soil[jsat][l] - theta[l]) * hs * Qsat / r->co[l];
				evap_layer[l] = rho*(soil[sy][jsat][l] - theta[l]) * hs * Qsat / r[l];
				
				//A += ( (soil[jsat][l] - theta[l]) / (soil[jsat][1] - theta[1]) ) * (rv/r->co[l]) * hs * Qsat;
				A += ( (soil[sy][jsat][l] - theta[l]) / (soil[sy][jsat][1] - theta[1]) ) * (rv/r[l]) * hs * Qsat;
				//B += ( (soil[jsat][l] - theta[l]) / (soil[jsat][1] - theta[1]) ) * (rv/r->co[l]);
				B += ( (soil[sy][jsat][l] - theta[l]) / (soil[sy][jsat][1] - theta[1]) ) * (rv/r[l]);
				
			}
			
			Qs = ( Qa + A ) / ( 1. + B );
			
			for(l=1;l<n;l++){
			//	evap_layer->co[l] -= rho*(soil[jsat][l] - theta[l]) * Qs / r->co[l];
				evap_layer[l] -= rho*(soil[sy][jsat][l] - theta[l]) * Qs / r[l];
			}
			
			//free_doublevector(r);
			
			//calculates evaporation from the surface
			F = ( soil[sy][jsat][1] ) / ( soil[sy][jsat][1] - soil[sy][jres][1] );
		//	evap_layer->co[1] += rho*( theta[1] - soil[jres][1] ) * F * ( Qgsat - Qa ) / rv;
			evap_layer[1] += rho*( theta[1] - soil[sy][jres][1] ) * F * ( Qgsat - Qa ) / rv;
															
			*beta = ( soil[sy][jsat][1] - theta[1] ) + ( theta[1] - soil[sy][jres][1] ) * F  - ( soil[sy][jsat][1] - theta[1] ) / ( 1. + B );
			*alpha = ( ( theta[1] - soil[sy][jres][1] ) * F + ( soil[sy][jsat][1] - theta[1] ) * A / ( Qgsat * (1. + B ) ) ) / (*beta);
		}
	}
}


