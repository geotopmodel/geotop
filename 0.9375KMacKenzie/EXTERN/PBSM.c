#include "constant.h"
#include "struct.geotop.09375.h"
#include "snow.09375.h"
#include "PBSM.h"

extern T_INIT *UV;
extern double NoV;

#define fluidtodyn 0.8
#define incr_V_factor 1


/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void BlowingSnow(long r, long c, float FetchDist, SNOW *snow, double T, double RH, double Uten, double StubHt, double *Transport, double *Sublimation, double *LatentH){

	float EffectiveStubHt, Ustar_main, Znod, Lambda_main, Beta_main, Ustn_main, Uten_Prob;
	float P, Ut;
	float BoundaryHt = 5.0;
	int Snow_Age;
	float Trans,Subl,Lat;
	double Wliq,Wice;

	EffectiveStubHt = (float)StubHt;
	if(EffectiveStubHt<1.E-6) EffectiveStubHt=1.E-6;

	Ustar_main = (float)(0.02264*pow(Uten,1.295)); //Eq. 6.2 rev.,  Ustar over fallow

	if (EffectiveStubHt > 0.01){
		Znod = (pow(Ustar_main,2.)/163.3) + 0.5*320*EffectiveStubHt*0.003; //Eq. 29, Snowcover Book
		Lambda_main = 320.0*3.0/1000.0*EffectiveStubHt;  //{Raupach Eq. 1}
		Beta_main = 170.0;
		Ustn_main  = Ustar_main*sqrt((Beta_main*Lambda_main)/(1+Beta_main*Lambda_main));
		Uten_Prob = (log(10.0/Znod))/ka *(Ustar_main - Ustn_main);
	}else{
		Uten_Prob = (float)Uten;
	}

	Snow_Age = (int)(snow->dimens_age->co[r][c]/3600.0);	//snowage in hours
	if(Snow_Age<1) Snow_Age=1;

	if(snow->lnum->co[r][c]>0){
		Wliq=snow->w_liq->co[snow->lnum->co[r][c]][r][c];
		Wice=snow->w_ice->co[snow->lnum->co[r][c]][r][c];
	}else{
		Wliq=0.0;
		Wice=0.0;
	}

	Probability_Threshold(Wliq, Wice, Snow_Age, snow->Psnow->co[r][c], T, Uten_Prob, &P, &Ut);

	Ut *= fluidtodyn;	//{Fluid threshold to Dynamical threshold}

	Pbsm (EffectiveStubHt, BoundaryHt, FetchDist, Ut, Uten, (float)T, (float)RH, &Trans, &Subl, &Lat);

//  {Probability, P, used to calculate Transport rate, Sublimation, and Latent Heat Flux}
	Trans*=P;		//{kg/m-width/s}
	Subl*=P;		//{kg/m2/s}
	Lat*=P;			//{J/m2/s}

	*Transport = (double)Trans;
	*Sublimation = (double)Subl;
	*LatentH = (double)Lat;
}


//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//Prarie Blowing Snow Model - Pomeroy et al. (1993)

void Pbsm (float Meht, float Fht, float Fetch, float Uthr, float Uten, float Temp, float Rel_H, float *Trans, float *Subl, float *LatentH){

//     Modified Calculations for Mean Particle Mass in this version
//     program to calculate blowing snow horizontal flux, sublimation rate
//     and latent heat flux due to snow sublimation for a variety of windspeeds,
//     boundary layers and surface conditions.

//     All variable and constants entered into the programme are in SI and
//     use Canadian Atmospheric Environement Service Meteorological data
//     format.  Snow transport is in kg per square meter per second
//     from the surface to 5 metres height.  Sublimation is totaled to the top
//     of the boundary layer for diffusion, based on the meteorological
//     Fetch and is expressed in millimeters of blowing snow lost over
//     a square meter of snow surface per half hour

  float   A,      Alpha,  B,      Bd,     Bound,  C,
  Diff,   DmDt,   Es,     H,
  Htran,  Hsalt,  Inc,    Lamb,   Lambda, Lb,
  Mpm,    Mpr,    Nh,     Nsalt,
  Nz,     Nuss,   Omega,  TQsalt, FQsalt,
  FQsum,  TQsum,  Qz,     RauTerm,
  Reyn,   SBsalt, Sbz,    SBsum,
  SigmaZ, Sigma2, SvDens, Usthr,  Ustar,
  UstarZ, Uz,     Vs,     Vsalt,  Sigma,
  Vsusp,  Z,      Zr,     Zstb;

  short a, end=0;

//     define constants for equations

 const float Qstar = 120;    //{Solar Radiation Input}
 const float M = 18.01;      //{molecular weight of water (kg/kmole)}
 const float R = 8313;       //{universal gas constant (J/(kmole K))}
 const float LATH = 2.838E6; //{latent heat of sublimation (J/kg) List 1949}
 const float DICE=900;       //{density of ice, assumed equal to blowing snow part (kg/m3)}
 const float ZD=0.3;         //{height of boundary-layer at xd (m) Takeuchi (1980)}
 const float XD=300;         //{Fetch to develop b-l to ZD (m)}
 const float rho = 1.2;      //{kg/m2}
 const float Beta = 170;     //{Cr/Cs = 170}
 const float C1 = 2.8;       //{2.3}
 const float C2 = 1.6;
 const float C3 = 4.2;       //{3.25} {e = 1/(C3*Ustar)}
 const float PI =3.14159;

//Compute stubble coefficients

   Zstb = 0.0048*Meht*100.0;      //{Lettau, used for susp Z0''}
   Lambda = 320.0*3.0/1000.0*Meht;    //{Raupach Eq. 1}
   Sigma = (PI*0.003)/(4.0*Meht); //{Raupach Eq. 4}

// Initialization

   TQsalt = 0.0; //{Total saltation flux}
   FQsalt = 0.0; //{Fence ht saltation flux}
   TQsum = 0.0;  //{Total Suspension}
   FQsum = 0.0;  //{Fence}
   SBsalt = 0.0;
   SBsum = 0.0;

   Temp  = Temp + tk; //Convert to Deg. K}

// Check for data errors    Fluxes set to zero for the hour}

   Lamb = (0.00063*Temp + 0.0673)/1.0; //{therm. cond. of atm. (J/(msK))}

   Diff = 2.06E-5*pow(Temp/273.0, 1.75); //{diffus. of w.vap. atmos. (m2/s}

   B = LATH * M/(R * Temp) - 1.0;


// find undersaturation of w. vapour at 2 metres}

   Es = 611.15 * exp(22.452*(Temp - 273.0)/Temp);  //{sat pressure}
   SvDens = (Es*M)/(R*Temp);                     //{sat density}
   Sigma2 = Rel_H - 1.0;           //{undersaturation at 2 m}

   if(Uten > Uthr) {

// define saltation parameters and calculate saltation
//  rate using 10/1987 MODEL OF BLOWING SNOW EQUATIONS}

		Usthr = 0.03697*Uthr;              //{Eq. 6.3}
		Ustar = 0.02264*pow(Uten, 1.295); //{Eq. 6.2 rev}
		RauTerm  = 1.0/((1.0-Sigma*Lambda)*(1.0+Beta*Lambda)); //{Raupach}

		Hsalt = C2/(2.0*g)*pow(Ustar,2.0);    //{Eq. 4.13}
		Nsalt = 2.0*rho/(C2*C3*Ustar)*(RauTerm - pow(Usthr,2.0)/pow(Ustar,2.0)); //{Eq. 4.14 updated}

		if(Nsalt > 0){

			TQsalt = C1*rho*Usthr/(g*C3*Ustar)*(pow(Ustar,2.0)*RauTerm - pow(Usthr,2.0)); //{Eq. 4.20}
			FQsalt = TQsalt;

// {Some saltation goes over the fence}
			if(Hsalt > Fht) {FQsalt = TQsalt*Fht/Hsalt;}

// {calculate sublimation rate in the saltation layer}

			Mpr= 100E-6;
			Htran = 0.9 * PI * pow(Mpr,2.0) * Qstar;
			Alpha = 5.0;

			SigmaZ = Sigma2 * (1.019 + 0.027 * log(Hsalt)); //{ Eq. 6.20, Revised in May. 1997}
			if(SigmaZ > -0.01) {SigmaZ = -0.01;}
			Vsalt = 0.6325 * Ustar + 2.3 * Usthr; //{Eq. 6.25}

			Reyn = (2.0 * Mpr * Vsalt)/1.88E-5;     //{Eq. 6.22}
			Nuss = 1.79 + 0.606 * sqrt(Reyn);     //{Eq. 6.21}
			A = Lamb * Temp * Nuss;
			C = 1.0/(Diff * SvDens * Nuss);
			DmDt = ((2.0 * PI * Mpr * SigmaZ) - (Htran * B/A))/((LATH * B/A) + C);
			//{Eq. 6.16} {Gamma Dist. Corr.}
			Mpm = 4.0/3.0 * PI * DICE * Mpr*pow(Mpr,2.0) *(1.0 + 3.0/Alpha + 2.0/pow(Alpha,2.0));

			Vs = DmDt/Mpm;              //{Sublimation rate coefficient Eq. 6.13}

			SBsalt = Vs * Nsalt * Hsalt;  //{Eq. 6.11}

// calculate mass flux in the suspended layers and the sublimation
//   rate for layers of height Inc from height r to b}

			Zr = 0.05628 * Ustar;  //{Eq. 5.27}
			Alpha = 15.0;
			Inc = 0.01;

// Loop to find the first suspended drift density level, r
//   from the reference level Zr
//   To preserve continuity with saltation the first suspended
//   level drift density is less than or equal to Nsalt.}

			FQsum = 0;
			TQsum = 0;
			SBsum = 0;

			Z = Zr + Inc;

			do{
				if(Z<=0.15){
					Nz = 0.8 * exp(-1.55*(pow(0.05628*Ustar, -0.544) - pow(Z, -0.544))); // {Eq. 5.26, Revised in Dec. 1995}
					if(Nz>Nsalt) Z = Z + Inc;
				}
			}while(Z<=0.15 && Nz>Nsalt);

			Lb = Z + Inc;
			Z = Lb;
			Inc = 0.001;

// find height of fully-developed boundary layer for turbulent
//   diffusion using a form of Pasquills plume dispersion eq.
//   iterate towards Bound}

			Bd = 1.0;
			Bound = ZD + (ka*ka * (Fetch - XD) * pow(log(Bd * 162.926/
				pow(Ustar,2.0)) * log(ZD * 162.926/pow(Ustar,2.0)), -0.5));     //{Eq. 6.9}

			while (abs(Bound - Bd) > 0.001) {
				Bd = Bound;
				Bound = ZD + (ka*ka * (Fetch - XD) * pow(log(Bd * 162.926/
					pow(Ustar,2.0)) * log(ZD * 162.926/pow(Ustar,2.0)), -0.5)); //{Eq. 6.9}
			} //{while}

// Loop to calculate the suspended mass flux up to 5 metres
//   and the total sublimation rate to the top of the boundary layer
//   at increments of 1 mm to 50cm & increments of 10 cm to  b}

			do{
				a=0;
				H = Z + Inc;
				do{
					if (H <= Bound) {
						Nh = 0.8 * exp(-1.55*(pow(0.05628*Ustar, -0.544) - pow(H, -0.544)));
						Nz = Nh;
						UstarZ = Ustar * pow(1.2/(1.2 + Nz), 0.5);         //{Eq. 5.17a}
						Uz = (UstarZ/ka) *log(H/((0.00613 *pow(Ustar,2.0)) + Zstb));//{Eq. 4.17r}
						if(Uz > 0) {
							Mpr = 4.6E-5 * pow(H, -0.258);                     //{Eq. 6.15}
							if(H >= 5.0) {Mpr = 30E-6;}

							Htran = 0.9 * PI * pow(Mpr,2.0) * Qstar;
							Alpha = 4.08 + 12.6 * H;                             //{Eq. 6.14}

							if(H >= 1.5) {Alpha = 25.0;}

							SigmaZ = Sigma2 * (1.019 + 0.027 * log(H)); //{ Eq. 6.20, Revised in May. 1997}
							if(SigmaZ > -0.01) {SigmaZ = -0.01;}
							Omega = 1.1E7 * pow(Mpr, 1.8);              //{Eq. 5.18}
							Vsusp = Omega + 0.0106 * pow(Uz, 1.36);
							Reyn = (2.0 * Mpr * Vsusp)/1.88E-5;           //{Eq. 6.22}

							Nuss = 1.79 + 0.606 * sqrt(Reyn);           //{Eq. 6.21}
							A = Lamb * Temp * Nuss;
							C = 1.0/(Diff * SvDens * Nuss);
							DmDt = ((2.0*PI * Mpr * SigmaZ) - (Htran*B/A))/((LATH*B/A) + C);
							Mpm = 1.333 * PI * DICE * pow(Mpr,3.0) * (1.0 + 3.0/Alpha + 2.0/pow(Alpha,2.0));  //{Eq. 6.16} {Gamma Dist. Corr.}

							Vs = DmDt/Mpm;                             //{Eq. 6.13}

							Sbz = Vs * Nz * Inc;                       //{mg}
							SBsum = SBsum + Sbz;                       //{Eq. 6.12}
							Qz = Nz * Uz * Inc;                        //{Eq. 5.4}
							if(H >= 5.0) {Qz = 0.0;}
							TQsum = TQsum + Qz;                        //{Eq. 5.5}
							if(H <= Fht) {FQsum = FQsum + Qz;}         //{Eq. 5.5}

							if(Nz >= 1E-5) {
								if(((H-Inc) >= 0.5) && (H < 0.6)) {
									Inc = 0.1;
									Z = 0.5;
									a=1;   //{re start the loop}
								} //{if}
							}else{
								end=1;
							} //{if}

						}else {
							FQsalt = 0.0;        //{No wind at boundary; no flux}
							TQsalt = 0.0;
							FQsum = 0.0;
							TQsum = 0.0;
							SBsalt = 0.0;
							SBsum = 0.0;
							end=1;
						} //{if}
					}
					H = H + Inc;

				}while(H <= Bound && a==0 && end==0);	//{while}
			}while(a==1);
		}  //{if}
	} //{if}

	Sum (TQsalt, TQsum, SBsum, SBsalt, Trans, Subl, LatentH);

} //{PBSM procedure}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

  void Sum(float TQsalt, float TQsusp, float SBsum, float SBsalt, float *TransFlux, float *SubFlux, float *HeatFlux){

  float Qsubl2;
  float L = 2.838E6; //!latent heat of sublimation (J/kg) List 1949

// {total sublimation}

   Qsubl2 = (SBsum + SBsalt)*(-1.E+6); //{-mgmm2s to b MILLIGRAMS PER M2 S}
   if  ((SBsum + SBsalt) == 0.0) {Qsubl2 = 0.0;} //{- 0.000 to 0.000}

   *SubFlux = Qsubl2/(1.E6);        //{mg/m2/s to kg/m2/s}

   *TransFlux = (TQsalt + TQsusp); //{kg/m-width/s}
   *HeatFlux=  *SubFlux * L;          //{J/m2/s}

} //{sum procedure}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void Probability_Threshold(double Wliq, double Wice, int Snow_Age, double Psnow, double Temperature, float Uten, float *Probability, float *Threshold){

//{Probability of blowing snow occurrence and threshold wind speeds determined
//              by ambient air temperature and snow age}

	float Wind, Mean, Variance, c;

	if (Wice==0 && Psnow==0.0) {   //{no snow available}
        *Threshold = 9.43 + 0.18 * Temperature + 0.0033 * pow(Temperature,2.0); //{m/s}
        *Probability = 0.0;

	}else if (Wliq/Wice>0.01) {	//{or wet snow remains on the ground}
        Mean = 21.0;
        Variance = 7.0;
        Wind = 0.0;
        *Probability = 0.0;

        while ((Wind <= Uten) && (Uten >7.0)) {    //{loop to calculate P.} // {wind < 7 m/s too weak for wet snow transport}
            Wind = Wind + 0.1;
            c = (-pow(Wind - Mean,2.0))/(2.0*pow(Variance,2.0));
            *Probability = *Probability + (1.0/(Variance * 2.5055)) * (exp(c)) * 0.1;
        } //{while do}

        *Threshold = 9.9;     //{m/s}

	}else if(Psnow>0) {  // {with concurrent snowfall: new dry snow}
        Mean = 0.365 * Temperature + 0.00706 * pow(Temperature,2.0) + 0.91 * log(Snow_Age) + 11.0;
        Variance = 0.145 * Temperature + 0.00196 * pow(Temperature,2.0) + 4.23;
        Wind = 0.0;
        *Probability = 0.0;

        while ((Wind <= Uten) && (Uten >= 3.0)) {   // {Wind < 3 m/s too weak for dry snow transport}
           Wind = Wind + 0.1;
           c = (-pow(Wind - Mean,2.0))/(2.0*pow(Variance,2.0));
           *Probability = *Probability + (1.0/(Variance * 2.5055)) * (exp(c)) * 0.1;
        } //{while do}

        *Threshold = 9.43 + 0.18 * Temperature + 0.0033 * pow(Temperature,2.0); //{m/s}

   }else{  // {without concurrent snowfall: old dry snow}
        Mean = 0.365 * Temperature + 0.00706 * pow(Temperature,2.0) + 0.91 * log(Snow_Age) + 11.0;
        Variance = 0.145 * Temperature + 0.00196 * pow(Temperature,2.0) + 4.23;
        Wind = 0.0;
        *Probability = 0.0;

        while ((Wind <= Uten) && (Uten >= 3.0)) {   // {Wind < 3 m/s too weak for dry snow transport}
           Wind = Wind + 0.1;
           c = (-pow(Wind - Mean,2.0))/(2.0*pow(Variance,2.0));
           *Probability = *Probability + (1.0/(Variance * 2.5055)) * (exp(c)) * 0.1;
        } //{while do}

        *Threshold = 9.43 + 0.18 * Temperature + 0.0033 * pow(Temperature,2.0); //{m/s}
   } //{if}

} //{Probability_threshold procedure}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
void set_windtrans_snow(SNOW *snow, METEO *met, LAND *land, PAR *par, double t){

	float FetchDist=(float)par->snow_fetch;//between 300 and 1000 m
	double DW,SWE,D,rho,Wsub,ActToPot,DWl;
	double dx,dy,RH,T,U,StubHt;
	double LatH;
	long i,r,c,nr,nc,ns;
	long r0,c0;
	long num_change;
	short lu;
	double Qtrans=0.0,Qsubl=0.0,Dsum=0.0,SWEsum=0.0,Uav=0.0;//initialization
	long l;
	double q;

	double Wsubl_tot=0.0,Wtrans_tot=0.0;
	double W, a4, b4, c4, d4, e4, f4, CR;// h;

	a4=2.66E-3;c4=0.04;d4=0.0884;e4=0.046;f4=400;

	nr=snow->T->nrh;
	nc=snow->T->nch;

	dx=UV->U->co[1];
	dy=UV->U->co[2];

	//Wtrans [kg/m2] , Qtrans [kg m-1 s-1]
	initialize_doublematrix(snow->Wtrans, 0.0);

	//call PBSM
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(land->LC->co[r][c]!=NoV){

				//vegetation height above snow surface (stubble)
				lu=(short)land->LC->co[r][c];
				D=DEPTH(r,c,snow->lnum,snow->Dzl);
				//StubHt in [m], D and land->ty->co[lu][jz0thressoil] in [mm]
				if(D > land->ty->co[lu][jz0thressoil]){
					StubHt=0.0;
				}else{
					StubHt=land->ty->co[lu][jHveg] - 1.E-3*D;
				}

				//met variable
				RH=met->RHgrid->co[r][c];
				T=met->Tgrid->co[r][c];
				U=met->Vgrid->co[r][c];

				//find the equilibrium fluxes
				BlowingSnow(r, c, FetchDist, snow, T, RH, incr_V_factor*U, StubHt, &(snow->Qtrans->co[r][c]), &(snow->Qsub->co[r][c]), &LatH);

				//calculation only needed
				Qtrans += (fabs(snow->Qtrans->co[r][c]))/(double)par->total_pixel;
				Qsubl += snow->Qsub->co[r][c]/(double)par->total_pixel;
				Dsum += D/(double)par->total_pixel;
				Uav += U/(double)par->total_pixel;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					SWEsum += (snow->w_ice->co[l][r][c] + snow->w_liq->co[l][r][c])/(double)par->total_pixel;
				}
			}
		}
	}

	set_no_value(snow->Qtrans, land->LC, NoV);
	extend_topography(snow->Qtrans, NoV);


	//if there is snow and blowing snow at the same time
	if( (Qsubl>0 || Qtrans>0) && Dsum>0){

		//wind in direction west-east

		for(r=1;r<=nr;r++){

			for(c=1;c<=nc;c++){
				snow->Qtrans_x->co[r][c]=fabs(snow->Qtrans->co[r][c]*(-sin(met->Vdir->co[r][c]*Pi/180.)));
			}

			initialize_longvector(snow->change_dir_wind,0);
			num_change=0;
			c=1;

			num_change++;
			c0=c;
			snow->change_dir_wind->co[num_change]=c;

			//printf("R:%ld c0:%ld numchange:%ld chdir:%ld\n",r,c0,num_change,snow->change_dir_wind->co[num_change]);

			do{
				c=c0;
				do{
					c++;
				}while( (-sin(met->Vdir->co[r][c]*Pi/180.))*(-sin(met->Vdir->co[r][c0]*Pi/180.))>0 && c<nc );

				num_change++;
				c0=c;
				snow->change_dir_wind->co[num_change]=c;
				//printf("R:%ld c0:%ld numchange:%ld chdir:%ld\n",r,c0,num_change,snow->change_dir_wind->co[num_change]);

			}while(c0<nc);

			for(i=1;i<num_change;i++){
				if( (-sin(met->Vdir->co[r][snow->change_dir_wind->co[i]]*Pi/180.)) > 0 ){
					if(snow->change_dir_wind->co[i]!=1){
						snow->Qtrans_x->co[r][snow->change_dir_wind->co[i]]=0.0;
					}else{
						//snow->Qtrans_x->co[r][1]=0.0;
					}
					//printf("+...i:%ld min:%ld max:%ld\n",i,snow->change_dir_wind->co[i],snow->change_dir_wind->co[i+1]);
					for(c=snow->change_dir_wind->co[i]+1;c<=snow->change_dir_wind->co[i+1];c++){
						if(snow->change_dir_wind->co[i+1]==nc || (snow->change_dir_wind->co[i+1]!=nc && c<snow->change_dir_wind->co[i+1])){
							//increasing potential snow transport
							if( fabs(snow->Qtrans->co[r][c]*(-sin(met->Vdir->co[r][c]*Pi/180.))) >= fabs(snow->Qtrans->co[r][c-1]*(-sin(met->Vdir->co[r][c-1]*Pi/180.))) ){
								q=snow->Qtrans_x->co[r][c];
								snow->Qtrans_x->co[r][c]=snow->Qtrans_x->co[r][c-1]+(3./FetchDist)*dx*(snow->Qtrans_x->co[r][c]-snow->Qtrans_x->co[r][c-1]);
								if(q<snow->Qtrans_x->co[r][c]) printf("err1 q:%e Q:%e %ld %ld %ld\n",q,snow->Qtrans_x->co[r][c],l,r,c);
							//decreasing potential snow transport
							}else{
								snow->Qtrans_x->co[r][c]=Fmin(snow->Qtrans_x->co[r][c-1],snow->Qtrans_x->co[r][c]);
							}
							snow->Wtrans->co[r][c] += ( snow->Qtrans_x->co[r][c-1] - snow->Qtrans_x->co[r][c] )*par->Dt/dx;
						}
					}

				}else{
					if(snow->change_dir_wind->co[i+1]!=nc){
						snow->Qtrans_x->co[r][snow->change_dir_wind->co[i+1]-1]=0.0;
					}else{
						//snow->Qtrans_x->co[r][nc]=0.0;
					}
					//printf("-...i:%ld min:%ld max:%ld\n",i,snow->change_dir_wind->co[i],snow->change_dir_wind->co[i+1]);
					for(c=snow->change_dir_wind->co[i+1]-1;c>=snow->change_dir_wind->co[i];c--){
						if(snow->change_dir_wind->co[i+1]==nc || (snow->change_dir_wind->co[i+1]!=nc && c<snow->change_dir_wind->co[i+1]-1)){
							if( fabs(snow->Qtrans->co[r][c]*(-sin(met->Vdir->co[r][c]*Pi/180.))) >= fabs(snow->Qtrans->co[r][c+1]*(-sin(met->Vdir->co[r][c+1]*Pi/180.))) ){
								q=snow->Qtrans_x->co[r][c];
								snow->Qtrans_x->co[r][c]=snow->Qtrans_x->co[r][c+1]+(3./FetchDist)*dx*(snow->Qtrans_x->co[r][c]-snow->Qtrans_x->co[r][c+1]);
								if(q<snow->Qtrans_x->co[r][c]) printf("err2 q:%e Q:%e %ld %ld %ld\n",q,snow->Qtrans_x->co[r][c],l,r,c);
							}else{
								snow->Qtrans_x->co[r][c]=Fmin(snow->Qtrans_x->co[r][c+1],snow->Qtrans_x->co[r][c]);
							}
							snow->Wtrans->co[r][c] += ( snow->Qtrans_x->co[r][c+1] - snow->Qtrans_x->co[r][c] )*par->Dt/dx;
						}
					}
				}
			}
		}


		//wind in direction south-north
		for(c=1;c<=nc;c++){

			for(r=1;r<=nr;r++){
				snow->Qtrans_y->co[r][c]=fabs(snow->Qtrans->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.)));
			}

			initialize_longvector(snow->change_dir_wind,0);
			num_change=0;
			r=1;

			num_change++;
			r0=r;
			snow->change_dir_wind->co[num_change]=r;

			//printf("C:%ld r0:%ld numchange:%ld chdir:%ld\n",c,r0,num_change,snow->change_dir_wind->co[num_change]);

			do{
				r=r0;
				do{
					r++;
				}while( (-cos(met->Vdir->co[r][c]*Pi/180.))*(-cos(met->Vdir->co[r][c0]*Pi/180.))>0 && r<nr );

				num_change++;
				r0=r;
				snow->change_dir_wind->co[num_change]=r;

				//printf("C:%ld r0:%ld numchange:%ld chdir:%ld\n",c,r0,num_change,snow->change_dir_wind->co[num_change]);

			}while(r0<nr);

			for(i=1;i<num_change;i++){
				if( (-cos(met->Vdir->co[snow->change_dir_wind->co[i]][c]*Pi/180.)) < 0 ){
					if(snow->change_dir_wind->co[i]!=1){
						snow->Qtrans_y->co[snow->change_dir_wind->co[i]][c]=0.0;
					}else{
						//snow->Qtrans_y->co[1][c]=0.0;
					}
					//printf("+...i:%ld min:%ld max:%ld\n",i,snow->change_dir_wind->co[i],snow->change_dir_wind->co[i+1]);
					for(r=snow->change_dir_wind->co[i]+1;r<=snow->change_dir_wind->co[i+1];r++){
						if(snow->change_dir_wind->co[i+1]==nr || (snow->change_dir_wind->co[i+1]!=nr && r<snow->change_dir_wind->co[i+1])){
							if(fabs(snow->Qtrans->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.)))>=fabs(snow->Qtrans->co[r-1][c]*(-cos(met->Vdir->co[r-1][c]*Pi/180.)))){
								q=snow->Qtrans_y->co[r][c];
								snow->Qtrans_y->co[r][c]=snow->Qtrans_y->co[r-1][c]+(3./FetchDist)*dy*(snow->Qtrans_y->co[r][c]-snow->Qtrans_y->co[r-1][c]);
								if(q<snow->Qtrans_y->co[r][c]) printf("err3 q:%e Q:%e - %e %e - %e - i:%ld %ld %ld - %ld %ld %ld\n",q,snow->Qtrans_y->co[r][c],
									fabs(snow->Qtrans->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.))),
									fabs(snow->Qtrans->co[r-1][c]*(-cos(met->Vdir->co[r-1][c]*Pi/180.))),
									snow->Qtrans_y->co[r-1][c],
									i,snow->change_dir_wind->co[i],snow->change_dir_wind->co[i+1],
									l,r,c);
							}else{
								snow->Qtrans_y->co[r][c]=Fmin(snow->Qtrans_y->co[r-1][c],snow->Qtrans_y->co[r][c]);
							}
							snow->Wtrans->co[r][c] += ( snow->Qtrans_y->co[r-1][c] - snow->Qtrans_y->co[r][c] )*par->Dt/dy;
						}
					}
					//printf("end");

				}else{
					if(snow->change_dir_wind->co[i+1]!=nr){
						snow->Qtrans_y->co[snow->change_dir_wind->co[i+1]-1][c]=0.0;
					}else{
						//snow->Qtrans_y->co[nr][c]=0.0;
					}
					//printf("+...i:%ld min:%ld max:%ld\n",i,snow->change_dir_wind->co[i+1],snow->change_dir_wind->co[i]);
					for(r=snow->change_dir_wind->co[i+1]-1;r>=snow->change_dir_wind->co[i];r--){
						if(snow->change_dir_wind->co[i+1]==nr || (snow->change_dir_wind->co[i+1]!=nr && r<snow->change_dir_wind->co[i+1]-1)){
							if(fabs(snow->Qtrans->co[r][c]*(-cos(met->Vdir->co[r][c]*Pi/180.)))>=fabs(snow->Qtrans->co[r+1][c]*(-cos(met->Vdir->co[r+1][c]*Pi/180.)))){
								q=snow->Qtrans_y->co[r][c];
								snow->Qtrans_y->co[r][c]=snow->Qtrans_y->co[r+1][c]+(3./FetchDist)*dy*(snow->Qtrans_y->co[r][c]-snow->Qtrans_y->co[r+1][c]);
								if(q<snow->Qtrans_y->co[r][c]) printf("err4 q:%e Q:%e %ld %ld %ld\n",q,snow->Qtrans_y->co[r][c],l,r,c);
							}else{
								snow->Qtrans_y->co[r][c]=Fmin(snow->Qtrans_y->co[r+1][c],snow->Qtrans_y->co[r][c]);
							}
							snow->Wtrans->co[r][c] += ( snow->Qtrans_y->co[r+1][c] - snow->Qtrans_y->co[r][c] )*par->Dt/dy;
						}
					}
					//printf("end");
				}
			}
		}

		//printf("OK\n");

		//update snow depth
		for(r=1;r<=nr;r++){
			for(c=1;c<=nc;c++){
				if(land->LC->co[r][c]!=UV->V->co[2]){

					ns=snow->lnum->co[r][c];

					D=DEPTH(r,c,snow->lnum,snow->Dzl);
					SWE=get_SWE(r,c,snow->lnum,snow->w_ice,snow->w_liq);

					if(snow->Qtrans->co[r][c]>0){
						ActToPot=(pow(pow(snow->Qtrans_x->co[r][c],2.0)+pow(snow->Qtrans_y->co[r][c],2.0),0.5))/snow->Qtrans->co[r][c];
					}else{
						ActToPot=0.0;
					}

					if(ActToPot>1.0001) printf("Error set_windtrans_snow %f %f %f\n",ActToPot,
						(pow(pow(snow->Qtrans_x->co[r][c],2.0)+pow(snow->Qtrans_y->co[r][c],2.0),0.5)),
						snow->Qtrans->co[r][c] );

					Wsub=-par->Dt*snow->Qsub->co[r][c]*ActToPot;

					Wtrans_tot+=snow->Wtrans->co[r][c]/(double)par->total_pixel;
					Wsubl_tot+=Wsub/(double)par->total_pixel;

					DW=snow->Wtrans->co[r][c] + Wsub;

					//snow compaction because of wind transport (Jordan, 99)
					W=0.0;
					for(l=snow->lnum->co[r][c];l>=1;l--){
						rho=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c])/(1.0E-3*snow->Dzl->co[l][r][c]);
						b4=Fmin(1.0,exp(-e4*(rho-f4)));
						W+=(snow->w_liq->co[l][r][c]+snow->w_ice->co[l][r][c]);	//kg m-2
						CR=-a4*b4*(pow(pow(snow->Qtrans_x->co[r][c],2.0)+pow(snow->Qtrans_y->co[r][c],2.0),0.5))*exp(-c4*snow->T->co[l][r][c]-d4*g*W);
						//snow->Dzl->co[l][r][c]*=(1.0 + CR*par->Dt);// was like this
						snow->Dzl->co[l][r][c] *= Fmax(0.1, 1.0 + CR*par->Dt);// new by stefano, 18/11/09
					}


					if(ns>0){	//snow on the soil

						if(DW<0){	//snow eroded

							i=ns;
							DWl=0.0;
							do{
								if(i<ns){

									if(snow->w_ice->co[i+1][r][c]<0){
										DW=snow->w_ice->co[i+1][r][c];
										DWl=snow->w_liq->co[i+1][r][c];
										snow->w_ice->co[i+1][r][c]=0.0;
										snow->w_liq->co[i+1][r][c]=0.0;
										snow->Dzl->co[i+1][r][c]=0.0;
										snow->lnum->co[r][c]-=1;
									}
								}

								snow->Dzl->co[i][r][c]*=(snow->w_ice->co[i][r][c]+DW)/snow->w_ice->co[i][r][c];

								/*if(snow->Dzl->co[i][r][c]>0){
									h = internal_energy(snow->w_ice->co[i][r][c], snow->w_liq->co[i][r][c], snow->T->co[i][r][c]);
								}*/

								snow->w_ice->co[i][r][c]+=DW;				//kg/m2
								snow->w_liq->co[i][r][c]+=DWl;				//kg/m2

								/*if(snow->Dzl->co[i][r][c]>0){
									from_internal_energy(r, c, h + internal_energy(DW, 0.0, snow->T->co[i][r][c]) + Lf*DWl, &(snow->w_ice->co[i][r][c]),
										&(snow->w_liq->co[i][r][c]), &(snow->T->co[i][r][c]));
								}*/

								DW=0.0;
								DWl=0.0;

								i--;

							}while(snow->w_ice->co[i+1][r][c]<0 && i>0);

							if(i==0 && snow->w_ice->co[i+1][r][c]<0){
								snow->w_ice->co[i+1][r][c]=0.0;				//kg/m2
								snow->Dzl->co[i+1][r][c]=0.0;	//mm
								snow->lnum->co[r][c]=0;
							}

						}else{	//snow drifted

							i = snow->lnum->co[r][c];

							//h = internal_energy(snow->w_ice->co[i][r][c], snow->w_liq->co[i][r][c], snow->T->co[i][r][c]);

							snow->w_ice->co[i][r][c]+=DW;

							/*from_internal_energy(r, c, h + internal_energy(DW, 0.0, snow->T->co[i][r][c]), &(snow->w_ice->co[i][r][c]),
								&(snow->w_liq->co[i][r][c]), &(snow->T->co[i][r][c]));*/

							snow->Dzl->co[snow->lnum->co[r][c]][r][c]+=1.0E+3*DW/rho_newlyfallensnow(met->Vgrid->co[r][c], snow->T->co[i][r][c], Tfreezing);

						}

					}else{	//snot not on the soil

						if(DW>0){

							snow->w_ice->co[1][r][c]+=DW;
							snow->Dzl->co[1][r][c]+=1.0E+3*DW/rho_newlyfallensnow(met->Vgrid->co[r][c], met->Tgrid->co[r][c], Tfreezing);
							snow->T->co[1][r][c]=T;
							if(snow->T->co[1][r][c]>Tfreezing) snow->T->co[1][r][c]=Tfreezing;
						}

					}

					snow_layer_combination(r, c, snow, T, par->snowlayer_inf, par->Dmin, par->Dmax, t);

				}
			}
		}
	}

	Dsum=0.0;
	SWEsum=0.0;
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(land->LC->co[r][c]!=NoV){
				//vegetation height above snow surface (stubble)
				D=DEPTH(r,c,snow->lnum,snow->Dzl);
				Dsum += D/(double)par->total_pixel;
				for(l=1;l<=snow->lnum->co[r][c];l++){
					SWEsum += (snow->w_ice->co[l][r][c] + snow->w_liq->co[l][r][c])/(double)par->total_pixel;
				}
			}
		}
	}
	printf("BLOWING SNOW: Dend:%f SWEend:%f U:%f Wsubl:%e Wtrans:%e DW:%e SWEin:%f\n",Dsum,SWEsum,Uav,Wsubl_tot,Wtrans_tot,Wsubl_tot+Wtrans_tot,SWEsum-Wsubl_tot-Wtrans_tot);

}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************

void print_windtrans_snow(long r, long c, SNOW *snow, PAR *par){

	long i;
	double ActToPot, Wsub, DW;

	if(snow->Qtrans->co[r][c]>0){
		ActToPot=(pow(pow(snow->Qtrans_x->co[r][c],2.0)+pow(snow->Qtrans_y->co[r][c],2.0),0.5))/snow->Qtrans->co[r][c];
	}else{
		ActToPot=0.0;
	}
	Wsub=-par->Dt*snow->Qsub->co[r][c]*ActToPot;
	DW=snow->Wtrans->co[r][c]+Wsub;

	if(par->output_snow>0){
		snow->Wtot->co[r][c]+=DW;
		//snow->Wtrans_cum->co[r][c]+=snow->Wtrans->co[r][c];
		//snow->Wsusp_cum->co[r][c]+=snow->Wsub->co[r][c];
		//snow->Wsubl_cum->co[r][c]+=snow->Wsub->co[r][c];
		//snow->Wsubgrid_cum->co[r][c]+=snow->Qtrans->co[r][c];
	}

	if(par->state_pixel==1){
		for(i=1;i<=par->chkpt->nrh;i++){
			if(r==par->rc->co[i][1] && c==par->rc->co[i][2]){
				snow->out_bs->co[1][i]+=snow->Wtrans->co[r][c];
				snow->out_bs->co[2][i]+=snow->Wtrans->co[r][c];
				snow->out_bs->co[3][i]+=Wsub;
				snow->out_bs->co[4][i]+=Wsub;
				/*snow->out_bs->co[5][i]+=snow->Qsub->co[r][c];
				snow->out_bs->co[6][i]+=snow->Qsub->co[r][c];
				snow->out_bs->co[7][i]+=snow->Qtrans->co[r][c];
				snow->out_bs->co[8][i]+=snow->Qtrans->co[r][c];*/
				snow->out_bs->co[9][i]+=DW;
				snow->out_bs->co[10][i]+=DW;
			}
		}
	}
}

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
void extend_topography(DOUBLEMATRIX *M, double novalue){

	long r,c,rr,cc;

	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(M->co[r][c]==novalue){
				find_the_nearest(r, c, novalue, M, &rr, &cc);
				M->co[r][c]=M->co[rr][cc];
			}
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

void find_the_nearest(long r, long c, double novalue, DOUBLEMATRIX *M, long *rr, long *cc){

	long i=0;
	short k;

	do{
		i++;
		k=0;
		if(k==0) k=no_novalue(r-i, c,   M, novalue, rr, cc);
		if(k==0) k=no_novalue(r-i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r  , c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c+i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c,   M, novalue, rr, cc);
		if(k==0) k=no_novalue(r+i, c-i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r,   c-i, M, novalue, rr, cc);
		if(k==0) k=no_novalue(r-i, c-i, M, novalue, rr, cc);
	}while(k==0);

}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/

short no_novalue(long r, long c, DOUBLEMATRIX *M, double novalue, long *rr, long *cc){

	short k=0;

	if((r>=1 && r<=M->nrh) && (c>=1 && c<=M->nch)){
		if(M->co[r][c]!=novalue){
			k=1;
			*rr=r;
			*cc=c;
		}
	}

	return(k);
}
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
void set_no_value(DOUBLEMATRIX *M, DOUBLEMATRIX *N, double undef){

	long r, c;

	for(r=1; r<=M->nrh; r++){
		for(c=1; c<=M->nch; c++){
			if(N->co[r][c]==undef) M->co[r][c]=undef;
		}
	}
}

/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
/******************************************************************************************************************************************/
