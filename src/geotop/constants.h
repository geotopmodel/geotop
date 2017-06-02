
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


//****************************************************
// Name of the program and no values 
//****************************************************

#define program_name "geotop.inpts"
#define logfile_name "geotop.log"

//****************************************************
// Keywords related
//****************************************************

#define max_charstring 200000
#define max_numvect 200000
#define num_par_number 415
#define num_par_char 367

//****************************************************
// Fixed Parameters
//****************************************************

#define KNe 0.0		//Euler method parameter for heat equation (0.5 = Crank Nicholson, 0 = Backward Euler)
#define LSAIthres 0.1 //Minimum LSAI 
#define z_evap 100. //soil depth responsable for soil evaporation [mm]
#define z_transp 10000. //soil depth responsable for canopy transpiration [mm]
#define min_tau_cloud 0.1
#define RelativeErrorRichards 1.E-10
#define max_cols_time_steps_file 100
#define PsiMin -1.E+10
#define thmin 0.1					//Newton method parameter 
#define thmax 0.5					//Newton method parameter 
#define Tol_h_mount 2.0 
#define Tol_h_flat 12.0
#define max_slope 89.999

//STANDARD LAPSE RATES: For the sign, remember that the Lapse Rate gives how a variable decrease with height 
#define LapseRateTair 6.5			//Lapse rate for Tair [C/m]
#define LapseRateTdew 2.5			//Lapse rate for Tdew [C/m]
#define LapseRatePrec 0.0			//Lapse rate for Precipitation [1/m]

//****************************************************
// Constants
//****************************************************

#define omega 0.261799388			/* Earth Rotation Velocity [rad/h] */
#define Isc 1367					/* Solar Constant  [W/m2] */
#define Pa0	1013.25					/* Mean atmospheric at sea level [mbar] */
#define rho_w 1000					/* density of water [kg/m3] */
#define rho_i 917					/* density of ice [kg/m3] */
#define Lf 333700.00				/* latent heat of fusion [J/kg] */
#define Ls 2834700.00				/* latent heat of sublimation [J/kg] */
#define g 9.81						/* gravity acceleration [m/s2] */
#define Pi 3.14159265358979			
#define tk 273.15					/* =0 Deg in Kelvin*/
#define k_liq 0.567					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_ice 2.290					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_air 0.023					/* thermal conductivity of air   [W m^-1 K^-1]*/
#define c_liq 4188.0				/* heat capacity of water		[J/(kg*K)]*/
#define c_ice 2117.0				/* heat capacity of ice		[J/(kg*K)]*/
#define c_can 2700.0				/* heat capacity of canopy [J/(kg*K)]*/
#define Tfreezing 0.0E1				/* freezing temperature [Celsius]*/
#define ka 0.41						/* Von Karman constant */
#define mu_l 0.001787				/* Dynamic viscosity of water at 0 degrees Celsius*/
#define wsn_vis 0.8					//snow on canopy: scattering parameters
#define wsn_nir 0.4
#define Bsnd_vis 0.5
#define Bsnd_nir 0.5
#define Bsnb_vis 0.5
#define Bsnb_nir 0.5
#define Cd 0.61						//discharge coefficient of the weir model
#define Rwv 461.495					//Specific gas constant for water vapor, 461.495 J/(kg·K)
#define D00 21.7					//molecular diffusivity of water vapor, 21.7 mm2/s
#define	secinday 86400.0			//seconds in one day

//****************************************************
//Meteo data
//****************************************************

#define iDate12	0					//Date12 : DDMMYYYYhhmm
#define iJDfrom0 iDate12+1			//Julian Day from year 0
#define iPrecInt iJDfrom0+1			//Precipitation
#define iPrec iPrecInt+1
#define iWs iPrec+1					//Total wind speed
#define iWdir iWs+1					//Wind direction
#define iWsx iWdir+1				//Wind speed from west, to east
#define iWsy iWsx+1					//Wind speed from south, to north
#define iRh iWsy+1					//Relative humidity
#define iT iRh+1					//Air temperature
#define iTdew iT+1					//Air dew temperature
#define iSW iTdew+1					//Global shortwave radiation
#define iSWb iSW+1					//Direct SW
#define iSWd iSWb+1					//Diffuse SW
#define itauC iSWd+1				//Cloud transmissivity in SWin
#define iC	 itauC+1				//Cloudiness factor
#define iLWi iC+1					//Incoming longwave
#define iSWn iLWi+1					//Net shortwave
#define iTs iSWn+1					//Surface Temperature
#define iTbottom iTs+1				//Bottom Temperature
#define nmet iTbottom+1	

//****************************************************
//soil data
//****************************************************

#define ilsDate12 0
#define ilsTa ilsDate12+1
#define ilsTdew	ilsTa+1
#define ilsPrec ilsTdew+1
#define nlstot ilsPrec+1

//****************************************************
//soil data
//****************************************************

#define jdz 1						//layer thickness [mm]
#define jpsi jdz+1					//initial psi [mm]
#define jT jpsi+1					//initial temperature [C]
#define jKn jT+1					//normal hydr. conductivity [mm/s]
#define jKl jKn+1					//lateral hydr. conductivity [mm/s]
#define jres jKl+1					//residual wat.cont.
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

//****************************************************
//land use data	
//****************************************************

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
#define ja_vis_dry jcd+1			//ground albedo in the visible spectrum dry
#define ja_nir_dry ja_vis_dry+1		//ground albedo in the near infrared spectrum dry
#define ja_vis_sat ja_nir_dry+1		//ground albedo in the visible spectrum saturated
#define ja_nir_sat ja_vis_sat+1		//ground albedo in the near infrared spectrum saturated
#define jemg ja_nir_sat+1			//soil emissivity
#define jcm jemg+1					//gauckler strickler (1/manning) coefficient
#define jN jcm+1					//number of roughness elements (vegetation) per unit surface for blowing snow
#define jdv	jN+1					//diameter of roughness elements (vegetation) [mm]
#define nlandprop jdv				//number of land use properties

//****************************************************
//vegetation files (the numbers must respect the same order as the block above)
//****************************************************

#define jdHveg 1					//vegetation height
#define jdz0thresveg jdHveg+1		//threshold on snow depth to change roughness length to snow covered values in vegetated area
#define jdz0thresveg2 jdz0thresveg+1
#define jdLSAI jdz0thresveg2+1		//LSAI
#define jdcf jdLSAI+1				//Canopy fraction
#define jddecay0 jdcf+1
#define jdexpveg jddecay0+1
#define jdroot jdexpveg+1
#define jdrs jdroot+1
#define jdvegprop jdrs

//****************************************************
//point output
//****************************************************

#define odate12 0
#define oJDfrom0 odate12+1
#define odaysfromstart oJDfrom0+1
#define operiod odaysfromstart+1
#define orun operiod+1
#define opoint orun+1
#define osnowover  opoint+1 //       prec_snow_atm;
#define orainover  osnowover+1  //        prec_rain_atm;
#define oprecsnow  orainover+1 //        prec_snow;
#define oprecrain  oprecsnow+1 //        (prec_rain_on_soil+prec_rain_on_snow);
#define orainonsnow oprecrain+1  //        prec_rain_on_snow;			
#define oV    orainonsnow+1   //          Vpoint/(double)n;
#define oVdir   oV+1  //     (met->Vdir->co[r][c])/(double)n;
#define oRH     oVdir+1  //    RHpoint/(double)n;
#define oPa     oRH+1  //   Ppoint/(double)n;
#define oTa     oPa+1 //    Tpoint/(double)n;
#define oTdew   oTa+1 //      Tdew/(double)n;
#define oTg     oTdew+1 //    Tg/(double)n; //Ts[C]
#define oTv     oTg+1 //    Tv/(double)n;
#define oTs     oTv+1 //    Ts/(double)n;
#define oEB     oTs+1 //    surfEB/(double)n; 
#define oG      oEB+1 //   G/(double)n; 
#define oSWin   oG+1 //      SWin/(double)n;
#define oSWb    oSWin+1 //     (SWbeam/(double)n);
#define oSWd    oSWb+1 //     (SWdiff/(double)n);
#define oLWin   oSWd+1 //      LWin/(double)n;	
#define ominLWin   oLWin+1 //      (epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
#define omaxLWin   ominLWin+1 //     (epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;
#define oSW     omaxLWin+1 //     SW/(double)n;
#define oLW     oSW+1 //     LW/(double)n;
#define oH      oLW+1 //    H/(double)n; //H[W/m^2]
#define oLE     oH+1 //     LE/(double)n; //ET[W/m^2]
#define ofc     oLE+1 //     fc/(double)n;
#define oLSAI    ofc+1 //      LSAI/(double)n;
#define oz0v    oLSAI+1 //      z0/(double)n;
#define od0v    oz0v+1 //      d0/(double)n;	
#define oEcan   od0v+1 //       (SWv+LWv-Hv-LEv)/(double)n;
#define oSWv    oEcan+1 //      SWv/(double)n;
#define oLWv    oSWv+1 //      LWv/(double)n;
#define oHv     oLWv+1 //     Hv/(double)n;
#define oLEv    oHv+1 //      LEv/(double)n;
#define oHg0    oLEv+1 //      Hg0/(double)n;
#define oLEg0   oHg0+1 //       Levap(Tg)*Eg0/(double)n;
#define oHg1    oLEg0+1 //      Hg1/(double)n;
#define oLEg1   oHg1+1 //       Levap(Tg)*Eg1/(double)n;
#define oevapsur  oLEg1+1 //        Er_soil*par->Dt;	//Eg[mm]
#define otrasp    oevapsur+1 //      Evt*(1000.0/rho_w)*par->Dt;	//Etc[mm]	
#define owcan_rain  otrasp+1 //        wat->wcan_rain->co[r][c]/(double)n;
#define owcan_snow  owcan_rain+1 //        wat->wcan_snow->co[r][c]/(double)n;	
#define oQv   owcan_snow+1 //       (Qv)/(double)n;
#define oQg   oQv+1 //       (Qg)/(double)n;
#define oQa   oQg+1 //       (Qa)/(double)n;
#define oQs   oQa+1 //       (Qs)/(double)n;
#define oLobuk  oQs+1 //        turbulence->co[2]/(double)n;
#define oLobukcan  oLobuk+1 //        (Locc)/(double)n;			
#define outop     oLobukcan+1 //     (u_top)/(double)n;
#define odecay     outop+1 //     (decay)/(double)n;
#define oSWup   odecay+1 //       SWupabove_v/(double)n;
#define oLWup     oSWup+1 //     LWup_above_v/(double)n;
#define oHup    oLWup+1 //      (H+fc*Hv)/(double)n;
#define oLEup   oHup+1 //       (LE+fc*LEv)/(double)n;
#define osnowdepth oLEup+1
#define oSWE osnowdepth+1
#define osnowdens oSWE+1
#define osnowT osnowdens+1
#define omrsnow   osnowT+1 //       Mr_snow*par->Dt;	//[mm]
#define osrsnow   omrsnow+1 //       Sr_snow*par->Dt;	//[mm]
#define oblowingsnowtrans osrsnow+1
#define oblowingsnowsubl oblowingsnowtrans+1
#define oglacdepth oblowingsnowsubl+1
#define oGWE oglacdepth+1
#define oglacdens oGWE+1
#define oglacT oglacdens+1
#define omrglac   oglacT+1 //       Mr_glac*par->Dt;	//[mm]
#define osrglac   omrglac+1 //       Sr_glac*par->Dt;	//[mm]
#define othawedup  osrglac+1 // thawed soil depth [mm]
#define othaweddw  othawedup+1 // thawed soil depth [mm]
#define owtableup  othaweddw+1 // water table depth [mm]
#define owtabledw  owtableup+1 // water table depth [mm]
#define otot owtabledw+1 // TOTAL NUMBER

//****************************************************
//BASIN OUTPUT
//****************************************************

#define oodate12 0
#define ooJDfrom0 oodate12+1
#define oodaysfromstart ooJDfrom0+1
#define ooperiod oodaysfromstart+1
#define oorun ooperiod+1
#define ooprecrain oorun+1										
#define ooprecsnow ooprecrain+1
#define oorainover ooprecsnow+1										
#define oosnowover oorainover+1
#define oopnet oosnowover+1
#define ooTa oopnet+1															
#define ooTg ooTa+1
#define ooTv ooTg+1			
#define ooevapsur ooTv+1															
#define ootrasp ooevapsur+1																	
#define ooLE ootrasp+1															
#define ooH ooLE+1															
#define ooSW ooH+1											
#define ooLW ooSW+1
#define ooLEv ooLW+1
#define ooHv ooLEv+1
#define ooSWv ooHv+1		
#define ooLWv ooSWv+1
#define ooSWin ooLWv+1
#define ooLWin ooSWin+1
#define oomasserror ooLWin+1
#define ootimestep oomasserror+1
#define ootot ootimestep+1 // TOTAL NUMBER

//****************************************************
//Files
//****************************************************

#define ftsteps 0			//file with time steps
#define fspar ftsteps+1				//soil parameters
#define fmet fspar+1				//meteo
#define fmetstlist fmet+1
#define fLRs fmetstlist+1					//lapse rates
#define fhormet fLRs+1					//horizon of meteo stations
#define fpointlist fhormet+1
#define fhorpoint fpointlist+1				//horizon of points for which the simulation is run 1D (point_sim==1)
#define fvegpar fhorpoint+1				//vegetation parameter
#define fqin fvegpar+1
#define fdem fqin+1					//digital elevation model (m)
#define flu fdem+1					//land use
#define fsoil flu+1					//soil type map
#define fdelay fsoil+1
#define fsky fdelay+1				//sky view factor
#define fslp fsky+1					//slope
#define fnet fslp+1				//channel network
#define fasp fnet+1					//aspect (0 north, then clockwise)
#define fcurv fasp+1				//curvature
#define fbed fcurv+1					//bedrock topography (m)
#define fwt0 fbed+1
#define fsn0 fwt0+1					//initial snow depth (mm)
#define fswe0 fsn0+1
#define fsnag0 fswe0+1				//initial snow age (days)
#define fgl0 fsnag0+1				//initial glacier depth (mm)
#define fQ fgl0+1					//(o.) output discharge file

#define fbas fQ+1					//o. basin variables
#define fbaswriteend fbas+1

#define fpoint fbaswriteend+1		//o. point variables
#define fpointwriteend fpoint+1

#define fTz fpointwriteend+1		//o. temperature profiles
#define fTzwriteend fTz+1

#define fTzav fTzwriteend+1
#define fTzavwriteend fTzav+1

#define fpsiz fTzavwriteend+1				//o. psi profiles
#define fpsizwriteend fpsiz+1

#define fpsiztot fpsizwriteend+1			//o. psi profiles
#define fpsiztotwriteend fpsiztot+1

#define fliqz fpsiztotwriteend+1			//o. water content profiles
#define fliqzwriteend fliqz+1

#define fliqzav fliqzwriteend+1			//o. water content profiles
#define fliqzavwriteend fliqzav+1

#define ficez fliqzavwriteend+1				//o. ice content profiles
#define ficezwriteend ficez+1

#define ficezav ficezwriteend+1				//o. ice content profiles
#define ficezavwriteend ficezav+1

#define fsatz ficezavwriteend+1

#define fsnTz fsatz+1				//o. snow data
#define fsnlz fsnTz+1		
#define fsniz fsnlz+1		
#define fsndz fsniz+1		

#define fsnTzwriteend fsndz+1				//o. snow data
#define fsnlzwriteend fsnTzwriteend+1		
#define fsnizwriteend fsnlzwriteend+1		
#define fsndzwriteend fsnizwriteend+1		

#define fglz fsndzwriteend+1					//o. glacier data
#define fglzwriteend fglz+1					//o. glacier data

#define fSCA fglzwriteend+1					//file giving the fraction of snow free areas and corresponding properties 

#define fTrun fSCA+1
#define fwrun fTrun+1
#define fdUrun fwrun+1	
#define fSWErun fdUrun+1
#define fTmaxrun fSWErun+1
#define fTminrun fTmaxrun+1
#define fwmaxrun fTminrun+1
#define fwminrun fwmaxrun+1

#define fT fwminrun+1					//o. temperature maps
#define fTsup fT+1					//o. temperature maps
#define fTav fTsup+1
#define fTavsup fTav+1
#define fliq fTavsup+1				//o. water content maps
#define fliqsup fliq+1				//o. water content in the soil at the surface
#define fliqav fliqsup+1
#define fice fliqav+1				//o. ice content maps
#define ficesup fice+1				//o. ice content maps
#define ficeav ficesup+1
#define fhsupland ficeav+1				//o. water over the surface (mm) maps
#define fhsupch fhsupland+1				//o. water over the surface (mm) maps

#define fradnet fhsupch+1				//o. radiation maps
#define fradLWin fradnet+1
#define fradLW fradLWin+1
#define fradSW fradLW+1
#define fradSWin fradSW+1
#define fradSWinbeam fradSWin+1
#define fshadow fradSWinbeam+1

#define fG fshadow+1				//o. surface heat flux maps
#define fH fG+1						//o. sensible heat flux maps
#define fLE fH+1					//o. latent heat flux maps
#define fTs fLE+1					//o. surface temperature maps
#define fprec fTs+1					//o. precipitation maps
#define fcint fprec+1				//o. precipitation intercepted by canopy maps
#define fpsiliq fcint+1				//o. psi maps
#define fpsitot fpsiliq+1
#define fsnowdepth fpsitot+1			//o. snow maps
#define fglacdepth fsnowdepth+1		//o. glacier maps
#define fsnowmelt fglacdepth+1		//o. snow melted maps
#define fsnowsubl fsnowmelt+1		//o. snow sublimated maps
#define fglacmelt fsnowsubl+1		//o. glacier ice melted maps
#define fglacsubl fglacmelt+1		//o. glacier ice sublimated maps
#define fTa fglacsubl+1				//o. air temperature maps
#define fwspd fTa+1					//o. wind speed maps
#define fwdir fwspd+1				//o. wind direction maps
#define frh fwdir+1					//o. relative humidity maps
#define fswe frh+1					//o. snow density maps
#define fgwe fswe+1					//o. glacier ice density maps
#define fsndur fgwe+1				//o. snow duration maps (hrs)
#define fthawed_up fsndur+1				//thawed soil depth
#define fthawed_dw fthawed_up+1
#define fwtable_up fthawed_dw+1
#define fwtable_dw fwtable_up+1			//water table depth
#define fpnet fwtable_dw+1
#define fevap fpnet+1

#define pG fevap+1
#define pH pG+1
#define pLE pH+1
#define pHg pLE+1					//specific day map plots(p.) sensible heat flux
#define pLEg pHg+1					//p. latent heat flux 
#define pHv pLEg+1
#define pLEv pHv+1
#define pSWin pLEv+1				//p. incoming shortwave radiation
#define pSWg pSWin+1				//p. outgoing shortwave radiation
#define pSWv pSWg+1
#define pLWin pSWv+1				//p. incoming longwave radiation
#define pLWg pLWin+1				//p. outgoing longwave radiation
#define pLWv pLWg+1
#define pTs pLWv+1					//p. surface temperature
#define pTg pTs+1
#define pTv pTg+1
#define pTa pTv+1					//p. air temperature
#define pVspd pTa+1					//p. wind speed
#define pVdir pVspd+1				//p. wind direction
#define pRH pVdir+1					//p. relative humidity
#define pD pRH+1					//p. snow depth
#define pth pD+1					//p. water content of the most superficial layer


#define rpsi pth+1					//recover file (f.) psi
#define riceg rpsi+1				//r. soil ice content
#define rTg riceg+1					//r. soil temperature
#define rDzs rTg+1					//r. snow layer thicknesses
#define rwls rDzs+1					//r. snow liquid water contents
#define rwis rwls+1					//r. snow ice contents
#define rTs rwis+1					//r. snow temperatures
#define rDzi rTs+1					//r. glacier layer thicknesses
#define rwli rDzi+1					//r. glacier liquid water contents
#define rwii rwli+1					//r. glacier ice contents
#define rTi rwii+1					//r. glacier temperatures
#define rns rTi+1						//r. number of snow layers
#define rni rns+1						//r. number of glacier layers
#define rsnag rni+1					//r. snow age
#define rwcrn rsnag+1				//r. water stored on canopy
#define rwcsn rwcrn+1
#define rTv rwcsn+1
#define rpsich rTv+1			
#define ricegch rpsich+1			
#define rTgch ricegch+1
#define rTrun rTgch+1
#define rwrun rTrun+1
#define rdUrun rwrun+1	
#define rSWErun rdUrun+1
#define rTmaxrun rSWErun+1
#define rTminrun rTmaxrun+1
#define rwmaxrun rTminrun+1
#define rwminrun rwmaxrun+1
#define rtime rwminrun+1
#define rsux rtime+1

#define nfiles rsux+1					//number of files

//****************************************************
//Points
//****************************************************
#define ptID 1
#define ptX ptID+1
#define ptY ptX+1
#define ptZ ptY+1
#define ptLC ptZ+1
#define ptSY ptLC+1
#define ptS ptSY+1
#define ptA ptS+1
#define ptSKY ptA+1
#define ptCNS ptSKY+1
#define ptCWE ptCNS+1
#define ptCNwSe ptCWE+1
#define ptCNeSw ptCNwSe+1
#define ptDrDEPTH ptCNeSw+1
#define ptHOR ptDrDEPTH+1
#define ptMAXSWE ptHOR+1
#define ptLAT ptMAXSWE+1
#define ptLON ptLAT+1
#define ptBED ptLON+1
#define ptTOT ptBED

