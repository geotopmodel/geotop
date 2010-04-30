
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


	
//Constants

/* Name of the program */

#define PROGRAM_NAME  "__geotop"


//****************************************************
//Constants to modify 
//****************************************************
#define LSAIthres 0.1
#define z_evap 40. //soil depth responsable for soil evaporation [mm]
#define z_transp 1000. //soil depth responsable for canopy transpiration [mm]

//****************************************************
//Constants to fixed 
//****************************************************
#define omega 0.261799388			/* velocita' di rotazione terrestre [rad/hr] */
#define Isc 1367					/* Costante solare [W/mq] */
#define Pa0	1013.25					/* Mean atmospheric at sea level [mbar] */
#define rho_w 1000					/* density of water [kg/mc] */
#define rho_i 917					/* density of ice [kg/mc] */
#define Lf 333700.00				/* heat of fusion [J/kg] */
#define g 9.81						/* gravity acceleration [m/s2] */
#define Pi 3.14159265358979			/* greek P*/
#define tk 273.15					/* =0 Deg in Kelvin*/
#define k_liq 0.567					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_ice 2.290					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_air 0.023					/* thermal conductivity of air   [W m^-1 K^-1]*/
#define c_liq 4188.0				/* heat capacity of water		[J/(kg*K)]*/
#define c_ice 2117.0				/* heat capacity of ice		[J/(kg*K)]*/
#define c_can 2700.0				/* heat capacity of canopy [J/(kg*K)]*/
#define KNe 0.5						/* Krank-Nicholson parameter for energy balance*/
#define Tfreezing 0.0E1				/* freezing temperature [Celsius]*/
#define ka 0.41						/* Von Karman constant*/
#define mu_l 0.001787				/* Dynamic viscosity of water at 0 degrees Celsius*/
#define Asurr 0.0
#define wsn_vis 0.8					//snow on canopy: scattering parameters
#define wsn_nir 0.4
#define Bsnd_vis 0.5
#define Bsnd_nir 0.5
#define Bsnb_vis 0.5
#define Bsnb_nir 0.5
#define PsiMin -1.E+10
#define Cd 0.61						//discharge coefficient of the weir model
#define thmin 0.1					//Newton method parameter 
#define thmax 0.5					//Newton method parameter 
#define Rwv 461.495					//Specific gas constant for water vapor, 461.495 J/(kg·K)
#define D00 21.7					//molecular diffusivity of water vapor, 21.7 mm2/s
#define	secinday 86400.0			//seconds in one day
#define RelativeErrorRichards 1.E-10
#define max_time_reduction_time 7

	
//Meteo data
#define iPt 0						/*Precipitation*/
#define iWs iPt+1					/*Wind speed*/
#define iWd iWs+1					/*Wind direction*/
#define iRh iWd+1					/*Relative humidity*/
#define iT iRh+1					/*Air temperature*/
#define iPs iT+1					/*Air Pressure*/
#define iSW iPs+1					/*global shortwave radiation*/
#define iSWb iSW+1					/*direct SW*/
#define iSWd iSWb+1					/*diffuse SW*/
#define itauC iSWd+1
#define iC	 itauC+1				/*Cloudiness*/
#define iLWi iC+1					/*incoming longwave*/
#define iSWn iLWi+1					/*net shortwave*/
#define nmet iSWn+1

//soil data
#define jdz 1						//layer thickness [mm]
#define jpsi jdz+1					//initial psi [mm]
#define jT jpsi+1					//initial temperature [C]
#define jKh jT+1					//lateral hydr. conductivity [mm/s]
#define jKv jKh+1					//vertical hydr. conductivity [mm/s]
#define jres jKv+1					//residual wat.cont.
#define jwp jres+1					//wilting point water cont.
#define jfc jwp+1					//field capacity water cont. 
#define jsat jfc+1					//porosity
#define ja jsat+1					//alpha[mm^-1]
#define jns ja+1					//n
#define jv jns+1					//v
#define jkt jv+1					//thermal conductivity
#define jct jkt+1					//thermal capacity
#define jss jct+1					//soil specific storativity
#define nsoilprop jss				//number of soil properties considered

//land use data	
#define jz0 1						//roughness length for soil
#define jz0thressoil jz0+1			//threshold on snow depth to change roughness length to snow covered values in soil area
#define jHveg jz0thressoil+1		//vegetation height
#define jz0thresveg jHveg+1			//threshold on snow depth to change roughness length to snow covered values in vegetated area
#define jz0thresveg2 jz0thresveg+1
#define jLSAI jz0thresveg2+1		//LSAI
#define jcf jLSAI+1					//Canopy fraction
#define jdecay0 jcf+1
#define jexpveg jdecay0+1
#define jroot jexpveg+1				//root depth [mm]
#define jrs jroot+1					//canopy transpiration coefficient
#define jvR_vis jrs+1				//vegetation in the visible spectrum
#define jvR_nir jvR_vis+1			//vegetation in the near infrared spectrum
#define jvT_vis jvR_nir+1			//vegetation in the visible spectrum
#define jvT_nir jvT_vis+1			//vegetation in the near infrared spectrum
#define jvCh jvT_nir+1				//departure of leaf angles from a random distribution (1 horizontal, 0 random, -1 vertical)
#define jcd jvCh+1					//surface density of canopy [kg/(m2*LSAI)]
#define ja_vis jcd+1				//ground albedo in the visible spectrum
#define ja_nir ja_vis+1				//ground albedo in the near infrared spectrum
#define jemg ja_nir+1				//soil emissivity
#define jcm jemg+1					//gauckler strickler (1/manning) coefficient
#define nlandprop jcm				//number of land use properties


//vegetation files (the numbers must respect the same order as the block above)
#define jdHveg 1					//vegetation height
#define jdz0thresveg jdHveg+1		//threshold on snow depth to change roughness length to snow covered values in vegetated area
#define jdz0thresveg2 jdz0thresveg+1
#define jdLSAI jdz0thresveg2+1		//LSAI
#define jdcf jdLSAI+1				//Canopy fraction
#define jddecay0 jdcf+1
#define jdexpveg jddecay0+1
#define jdroot jdexpveg+1
#define	jdrs jdroot+1
#define jdvegprop jdrs


//point output
#define osnowover  1 //       prec_snow_atm;
#define orainover  2 //        prec_rain_atm;
#define oprecsnow  3 //        prec_snow;
#define oprecrain  4 //        (prec_rain_on_soil+prec_rain_on_snow);
#define orainonsnow 5  //        prec_rain_on_snow;			
#define oV    6   //          Vpoint/(double)n;
#define oVdir   7  //     (met->Vdir->co[r][c])/(double)n;
#define oRH     8  //    RHpoint/(double)n;
#define oPa     9  //   Ppoint/(double)n;
#define oTa     10 //    Tpoint/(double)n;
#define oTdew   11 //      Tdew/(double)n;
#define oTg     12 //    Tg/(double)n; //Ts[C]
#define oTv     13 //    Tv/(double)n;
#define oTs     14 //    Ts/(double)n;
#define oEB     15 //    surfEB/(double)n; 
#define oG      16 //   G/(double)n; 
#define oSWin   17 //      SWin/(double)n;
#define oSWb    18 //     (SWbeam/(double)n);
#define oSWd    19 //     (SWdiff/(double)n);
#define oLWin   20 //      LWin/(double)n;	
#define ominLWin   21 //      (epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
#define omaxLWin   22 //     (epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
#define oSW     23 //     SW/(double)n;
#define oLW     24 //     LW/(double)n;
#define oH      25 //    H/(double)n; //H[W/m^2]
#define oLE     26 //     LE/(double)n; //ET[W/m^2]
#define ofc     27 //     fc/(double)n;
#define oLSAI    28 //      LSAI/(double)n;
#define oz0v    29 //      z0/(double)n;
#define od0v    30 //      d0/(double)n;	
#define oEcan   31 //       (SWv+LWv-Hv-LEv)/(double)n;
#define oSWv    32 //      SWv/(double)n;
#define oLWv    33 //      LWv/(double)n;
#define oHv     34 //     Hv/(double)n;
#define oLEv    35 //      LEv/(double)n;
#define oHg0    36 //      Hg0/(double)n;
#define oLEg0   37 //       Levap(Tg)*Eg0/(double)n;
#define oHg1    38 //      Hg1/(double)n;
#define oLEg1   39 //       Levap(Tg)*Eg1/(double)n;
#define oevapsur  40 //        Er_soil*par->Dt;	//Eg[mm]
#define otrasp    41 //      Evt*(1000.0/rho_w)*par->Dt;	//Etc[mm]	
#define owcan_rain  42 //        wat->wcan_rain->co[r][c]/(double)n;
#define owcan_snow  43 //        wat->wcan_snow->co[r][c]/(double)n;	
#define oQv   44 //       (Qv)/(double)n;
#define oQg   45 //       (Qg)/(double)n;
#define oQa   46 //       (Qa)/(double)n;
#define oQs   47 //       (Qs)/(double)n;
#define oLobuk  48 //        turbulence->co[2]/(double)n;
#define oLobukcan  49 //        (Locc)/(double)n;			
#define outop     50 //     (u_top)/(double)n;
#define odecay     51 //     (decay)/(double)n;
#define oSWup   52 //       SWup_above_v/(double)n;
#define oLWup     53 //     LWup_above_v/(double)n;
#define oHup    54 //      (H+fc*Hv)/(double)n;
#define oLEup   55 //       (LE+fc*LEv)/(double)n;
#define omrsnow   56 //       Mr_snow*par->Dt;	//[mm]
#define oersnow   57 //       Er_snow*par->Dt;	//[mm]
#define osrsnow   58 //       Sr_snow*par->Dt;	//[mm]
#define omrglac   59 //       Mr_glac*par->Dt;	//[mm]
#define oerglac   60 //       Er_glac*par->Dt;	//[mm]
#define osrglac   61 //       Sr_glac*par->Dt;	//[mm]
#define othawed  62 // thawed soil depth [mm]
#define owtable  63 // water table depth [mm]
#define otot 63 // TOTAL NUMBER

//BASIN OUTPUT
#define ooprecrain 1										
#define ooprecsnow 2
#define oorainover 3										
#define oosnowover 4
#define ooTa 5															
#define ooTg 6
#define ooTv 7			
#define ooevapsur 8																
#define ootrasp 9																	
#define ooLE 10															
#define ooH 11															
#define ooSW 12											
#define ooLW 13
#define ooLEv 14
#define ooHv 15
#define ooSWv 16		
#define ooLWv 17
#define ooSWin 18
#define ooLWin 19
#define oomasserror 20
#define ootot 20 // TOTAL NUMBER
