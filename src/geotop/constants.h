#ifndef _GEOTOP_CONSTANTS_H
#define _GEOTOP_CONSTANTS_H


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
#define num_par_number 406
#define num_par_char 367

//****************************************************
// Fixed Parameters
//****************************************************
namespace GTConst {
    const double KNe = 0.0;  // Euler method parameter for heat equation (0.5 = Crank Nicholson, 0 = Backward Euler)
    const double LSAIthres = 0.1;  // Minimum LSAI
    const double z_evap = 100.;  // soil depth responsable for soil evaporation [mm]
    const double z_transp = 10000.;  // soil depth responsable for canopy transpiration [mm]
    const double min_tau_cloud = 0.1;
    const double RelativeErrorRichards = 1.E-10;
    const double max_cols_time_steps_file = 100;
    const double PsiMin = -1.E+10;
    const double thmin = 0.1;  // Newton method parameter
    const double thmax = 0.5;  // Newton method parameter
    const double Tol_h_mount = 2.0;
    const double Tol_h_flat = 12.0;
    const double max_slope = 89.999;

    // STANDARD LAPSE RATES: For the sign, remember that the Lapse Rate gives how a variable decrease with height
    const double LapseRateTair = 6.5;  // Lapse rate for Tair [C/m]
    const double LapseRateTdew = 2.5;  // Lapse rate for Tdew [C/m]
    const double LapseRatePrec = 0.0;  // Lapse rate for Precipitation [1/m]
//****************************************************
// Constants
//****************************************************
    const double omega = 0.261799388; /* Earth Rotation Velocity [rad/h] */
    const double Isc = 1367;          /* Solar Constant  [W/m2] */
    const double Pa0 = 1013.25;       /* Mean atmospheric at sea level [mbar] */
    const double rho_w = 1000;        /* density of water [kg/m3] */
    const double rho_i = 917;         /* density of ice [kg/m3] */
    const double Lf = 333700.00;      /* latent heat of fusion [J/kg] */
    const double Ls = 2834700.00;     /* latent heat of sublimation [J/kg] */
    const double g = 9.81;      /* gravity acceleration [m/s2] */
    const double Pi = 3.14159265358979;
    const double tk = 273.15;    /* =0 Deg in Kelvin*/
    const double k_liq = 0.567;  /* thermal conductivity of water [W m^-1 K^-1]*/
    const double k_ice = 2.290;  /* thermal conductivity of water [W m^-1 K^-1]*/
    const double k_air = 0.023;  /* thermal conductivity of air   [W m^-1 K^-1]*/
    const double c_liq = 4188.0; /* heat capacity of water    [J/(kg*K)]*/
    const double c_ice = 2117.0; /* heat capacity of ice    [J/(kg*K)]*/
    const double c_can = 2700.0; /* heat capacity of canopy [J/(kg*K)]*/
    const double Tfreezing = 0.0E1; /* freezing temperature [Celsius]*/
    const double ka = 0.41;         /* Von Karman constant */
    const double mu_l = 0.001787; /* Dynamic viscosity of water at 0 degrees Celsius*/
    const double wsn_vis = 0.8; /* snow on canopy: scattering parameters */
    const double wsn_nir = 0.4;
    const double Bsnd_vis = 0.5;
    const double Bsnd_nir = 0.5;
    const double Bsnb_vis = 0.5;
    const double Bsnb_nir = 0.5;
    const double Cd = 0.61; /* discharge coefficient of the weir model */
    const double Rwv = 461.495;  // Specific gas constant for water vapor, 461.495 J/(kg·K)
    const double D00 = 21.7;  // molecular diffusivity of water vapor, 21.7 mm2/s
    const double secinday = 86400.0;  // seconds in one day
} // end namespace GTConst

//****************************************************
//Meteo data
//****************************************************
const unsigned int iDate12 = 0;             /*Date12 : DDMMYYYYhhmm*/
const unsigned int iJDfrom0 = iDate12 + 1;  /*Julian Day from year 0*/
const unsigned int iPrecInt = iJDfrom0 + 1; /*Precipitation*/
const unsigned int iPrec = iPrecInt + 1;
const unsigned int iWs = iPrec + 1;  /*Total wind speed*/
const unsigned int iWdir = iWs + 1;  /*Wind direction*/
const unsigned int iWsx = iWdir + 1; /*Wind speed from west, to east*/
const unsigned int iWsy = iWsx + 1;  /*Wind speed from south, to north*/
const unsigned int iRh = iWsy + 1;   /*Relative humidity*/
const unsigned int iT = iRh + 1;     /*Air temperature*/
const unsigned int iTdew = iT + 1;   /*Air dew temperature*/
const unsigned int iSW = iTdew + 1;  /*Global shortwave radiation*/
const unsigned int iSWb = iSW + 1;   /*Direct SW*/
const unsigned int iSWd = iSWb + 1;  /*Diffuse SW*/
const unsigned int itauC = iSWd + 1; /*Cloud transmissivity in SWin*/
const unsigned int iC = itauC + 1;   /*Cloudiness factor*/
const unsigned int iLWi = iC + 1;    /*Incoming longwave*/
const unsigned int iSWn = iLWi + 1;  /*Net shortwave*/
const unsigned int iTs = iSWn + 1;   /*Surface Temperature*/
const unsigned int iTbottom = iTs + 1;   /*Bottom Temperature (NOT present in GEOtop 2.1)*/
const unsigned int nmet = iTbottom + 1;

//****************************************************
//soil data
//****************************************************
const unsigned int ilsDate12 = 0;
const unsigned int ilsTa = ilsDate12 + 1;
const unsigned int ilsTdew = ilsTa + 1;
const unsigned int ilsPrec = ilsTdew + 1;
const unsigned int nlstot = ilsPrec + 1;

//****************************************************
//soil data
//****************************************************
const unsigned int jdz = 1;         /*layer thickness [mm]*/
const unsigned int jpsi = jdz + 1;  /*initial psi [mm]*/
const unsigned int jT = jpsi + 1;   /*initial temperature [C]*/
const unsigned int jKn = jT + 1;    /*normal hydr. conductivity [mm/s]*/
const unsigned int jKl = jKn + 1;   /*lateral hydr. conductivity [mm/s]*/
const unsigned int jres = jKl + 1;  /*residual wat.cont.*/
const unsigned int jwp = jres + 1;  /*wilting point water cont.*/
const unsigned int jfc = jwp + 1;   /*field capacity water cont. */
const unsigned int jsat = jfc + 1;  /*porosity*/
const unsigned int ja = jsat + 1;   /*alpha[mm^-1]*/
const unsigned int jns = ja + 1;    /*n*/
const unsigned int jv = jns + 1;    /*v*/
const unsigned int jkt = jv + 1;    /*thermal conductivity*/
const unsigned int jct = jkt + 1;   /*thermal capacity*/
const unsigned int jss = jct + 1;   /*soil specific storativity*/
const unsigned int nsoilprop = jss; /*number of soil properties considered*/

//****************************************************
//land use data
//****************************************************
const unsigned int jz0 = 1;     /*roughness length for soil*/
const unsigned int jz0thressoil = jz0 + 1;   /*threshold on snow depth to change roughness length to snow covered values in soil area*/
const unsigned int jHveg = jz0thressoil + 1; /*vegetation height*/
const unsigned int jz0thresveg = jHveg + 1; /*threshold on snow depth to change roughness length to snow covered values in vegetated area*/
const unsigned int jz0thresveg2 = jz0thresveg + 1;
const unsigned int jLSAI = jz0thresveg2 + 1; /*LSAI*/
const unsigned int jcf = jLSAI + 1;          /*Canopy fraction*/
const unsigned int jdecay0 = jcf + 1;
const unsigned int jexpveg = jdecay0 + 1;
const unsigned int jroot = jexpveg + 1; /*root depth [mm]*/
const unsigned int jrs = jroot + 1;     /*canopy transpiration coefficient*/
const unsigned int jvR_vis = jrs + 1;   /*vegetation in the visible spectrum*/
const unsigned int jvR_nir = jvR_vis + 1; /*vegetation in the near infrared spectrum*/
const unsigned int jvT_vis = jvR_nir + 1; /*vegetation in the visible spectrum*/
const unsigned int jvT_nir = jvT_vis + 1; /*vegetation in the near infrared spectrum*/
const unsigned int jvCh = jvT_nir + 1; /*departure of leaf angles from a random distribution (1 horizontal, 0 random, -1 vertical)*/
const unsigned int jcd = jvCh + 1; /*surface density of canopy [kg/(m2*LSAI)]*/
const unsigned int ja_vis_dry = jcd + 1; /*ground albedo in the visible spectrum dry*/
const unsigned int ja_nir_dry = ja_vis_dry + 1; /*ground albedo in the near infrared spectrum dry*/
const unsigned int ja_vis_sat = ja_nir_dry + 1; /*ground albedo in the visible spectrum saturated*/
const unsigned int ja_nir_sat = ja_vis_sat + 1; /*ground albedo in the near infrared spectrum saturated*/
const unsigned int jemg = ja_nir_sat + 1; /*soil emissivity*/
const unsigned int jcm = jemg + 1; /*gauckler strickler (1/manning) coefficient*/
const unsigned int jN = jcm + 1; /*number of roughness elements (vegetation) per unit surface for blowing snow*/
const unsigned int jdv = jN + 1; /*diameter of roughness elements (vegetation) [mm]*/
const unsigned int nlandprop = jdv; /*number of land use properties*/

//****************************************************
//vegetation files (the numbers must respect the same order as the block above)
//****************************************************
const unsigned int jdHveg = 1; /*vegetation height*/
const unsigned int jdz0thresveg = jdHveg + 1; /*threshold on snow depth to change roughness length to snow covered values in vegetated area*/
const unsigned int jdz0thresveg2 = jdz0thresveg + 1;
const unsigned int jdLSAI = jdz0thresveg2 + 1; /*LSAI*/
const unsigned int jdcf = jdLSAI + 1;          /*Canopy fraction*/
const unsigned int jddecay0 = jdcf + 1;
const unsigned int jdexpveg = jddecay0 + 1;
const unsigned int jdroot = jdexpveg + 1;
const unsigned int jdrs = jdroot + 1;
const unsigned int jdvegprop = jdrs;

//****************************************************
//point output
//****************************************************
const unsigned int odate12 = 0;
const unsigned int oJDfrom0 = odate12 + 1;
const unsigned int odaysfromstart = oJDfrom0 + 1;
const unsigned int operiod = odaysfromstart + 1;
const unsigned int orun = operiod + 1;
const unsigned int opoint = orun + 1;
const unsigned int osnowover = opoint + 1;    /*prec_snow_atm;*/
const unsigned int orainover = osnowover + 1; /*prec_rain_atm;*/
const unsigned int oprecsnow = orainover + 1; /*prec_snow;*/
const unsigned int oprecrain = oprecsnow + 1; /*(prec_rain_on_soil+prec_rain_on_snow);*/
const unsigned int orainonsnow = oprecrain + 1; /*prec_rain_on_snow;*/
const unsigned int oV = orainonsnow + 1;        /*Vpoint/(double)n;*/
const unsigned int oVdir = oV + 1;   /*(met->Vdir->co[r][c])/(double)n;*/
const unsigned int oRH = oVdir + 1;  /*RHpoint/(double)n;*/
const unsigned int oPa = oRH + 1;    /*Ppoint/(double)n;*/
const unsigned int oTa = oPa + 1;    /*Tpoint/(double)n;*/
const unsigned int oTdew = oTa + 1;  /*Tdew/(double)n;*/
const unsigned int oTg = oTdew + 1;  /*Tg/(double)n; //Ts[C]*/
const unsigned int oTv = oTg + 1;    /*Tv/(double)n;*/
const unsigned int oTs = oTv + 1;    /*Ts/(double)n;*/
const unsigned int oEB = oTs + 1;    /*surfEB/(double)n;*/
const unsigned int oG = oEB + 1;     /*G/(double)n;*/
const unsigned int oSWin = oG + 1;   /*SWin/(double)n;*/
const unsigned int oSWb = oSWin + 1; /*(SWbeam/(double)n);*/
const unsigned int oSWd = oSWb + 1;  /*(SWdiff/(double)n);*/
const unsigned int oLWin = oSWd + 1; /*LWin/(double)n;*/
const unsigned int ominLWin = oLWin + 1; /*(epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;*/
const unsigned int omaxLWin = ominLWin + 1; /*(epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;*/
const unsigned int oSW = omaxLWin + 1;   /*SW/(double)n;*/
const unsigned int oLW = oSW + 1;        /*LW/(double)n;*/
const unsigned int oH = oLW + 1;         /*H/(double)n; //H[W/m^2]*/
const unsigned int oLE = oH + 1;         /*LE/(double)n; //ET[W/m^2]*/
const unsigned int ofc = oLE + 1;        /*fc/(double)n;*/
const unsigned int oLSAI = ofc + 1;      /*LSAI/(double)n;*/
const unsigned int oz0v = oLSAI + 1;     /*z0/(double)n;*/
const unsigned int od0v = oz0v + 1;      /*d0/(double)n;*/
const unsigned int oEcan = od0v + 1;     /*(SWv+LWv-Hv-LEv)/(double)n;*/
const unsigned int oSWv = oEcan + 1;     /*SWv/(double)n;*/
const unsigned int oLWv = oSWv + 1;      /*LWv/(double)n;*/
const unsigned int oHv = oLWv + 1;       /*Hv/(double)n;*/
const unsigned int oLEv = oHv + 1;       /*LEv/(double)n;*/
const unsigned int oHg0 = oLEv + 1;      /*Hg0/(double)n;*/
const unsigned int oLEg0 = oHg0 + 1;     /*Levap(Tg)*Eg0/(double)n;*/
const unsigned int oHg1 = oLEg0 + 1;     /*Hg1/(double)n;*/
const unsigned int oLEg1 = oHg1 + 1;     /*Levap(Tg)*Eg1/(double)n;*/
const unsigned int oevapsur = oLEg1 + 1; /*Er_soil*par->Dt; //Eg[mm]*/
const unsigned int otrasp = oevapsur + 1; /*Evt*(1000.0/rho_w)*par->Dt; //Etc[mm]*/
const unsigned int owcan_rain = otrasp + 1; /*wat->wcan_rain->co[r][c]/(double)n;*/
const unsigned int owcan_snow = owcan_rain + 1;                        /*wat->wcan_snow->co[r][c]/(double)n;*/
const unsigned int oQv = owcan_snow + 1; /*(Qv)/(double)n;*/
const unsigned int oQg = oQv + 1;        /*(Qg)/(double)n;*/
const unsigned int oQa = oQg + 1;        /*(Qa)/(double)n;*/
const unsigned int oQs = oQa + 1;        /*(Qs)/(double)n;*/
const unsigned int oLobuk = oQs + 1;     /*turbulence->co[2]/(double)n;*/
const unsigned int oLobukcan = oLobuk + 1; /*(Locc)/(double)n;*/
const unsigned int outop = oLobukcan + 1;  /*(u_top)/(double)n;*/
const unsigned int odecay = outop + 1;     /*(decay)/(double)n;*/
const unsigned int oSWup = odecay + 1;     /*SWupabove_v/(double)n;*/
const unsigned int oLWup = oSWup + 1;      /*LWup_above_v/(double)n;*/
const unsigned int oHup = oLWup + 1;       /*(H+fc*Hv)/(double)n;*/
const unsigned int oLEup = oHup + 1;       /*(LE+fc*LEv)/(double)n;*/
const unsigned int osnowdepth = oLEup + 1;
const unsigned int oSWE = osnowdepth + 1;
const unsigned int osnowdens = oSWE + 1;
const unsigned int osnowT = osnowdens + 1;
const unsigned int omrsnow = osnowT + 1;  /*Mr_snow*par->Dt;  //[mm]*/
const unsigned int osrsnow = omrsnow + 1; /*Sr_snow*par->Dt;  //[mm]*/
const unsigned int oblowingsnowtrans = osrsnow + 1;
const unsigned int oblowingsnowsubl = oblowingsnowtrans + 1;
const unsigned int oglacdepth = oblowingsnowsubl + 1;
const unsigned int oGWE = oglacdepth + 1;
const unsigned int oglacdens = oGWE + 1;
const unsigned int oglacT = oglacdens + 1;
const unsigned int omrglac = oglacT + 1;      /*Mr_glac*par->Dt;  //[mm]*/
const unsigned int osrglac = omrglac + 1;     /*Sr_glac*par->Dt;  //[mm]*/
const unsigned int othawedup = osrglac + 1;   /*thawed soil depth [mm]*/
const unsigned int othaweddw = othawedup + 1; /*thawed soil depth [mm]*/
const unsigned int owtableup = othaweddw + 1; /*water table depth [mm]*/
const unsigned int owtabledw = owtableup + 1; /*water table depth [mm]*/
const unsigned int otot = owtabledw + 1;      /*TOTAL NUMBER*/

//****************************************************
//BASIN OUTPUT
//****************************************************
const unsigned int oodate12 = 0;
const unsigned int ooJDfrom0 = oodate12 + 1;
const unsigned int oodaysfromstart = ooJDfrom0 + 1;
const unsigned int ooperiod = oodaysfromstart + 1;
const unsigned int oorun = ooperiod + 1;
const unsigned int ooprecrain = oorun + 1;
const unsigned int ooprecsnow = ooprecrain + 1;
const unsigned int oorainover = ooprecsnow + 1;
const unsigned int oosnowover = oorainover + 1;
const unsigned int oopnet = oosnowover + 1;
const unsigned int ooTa = oopnet + 1;
const unsigned int ooTg = ooTa + 1;
const unsigned int ooTv = ooTg + 1;
const unsigned int ooevapsur = ooTv + 1;
const unsigned int ootrasp = ooevapsur + 1;
const unsigned int ooLE = ootrasp + 1;
const unsigned int ooH = ooLE + 1;
const unsigned int ooSW = ooH + 1;
const unsigned int ooLW = ooSW + 1;
const unsigned int ooLEv = ooLW + 1;
const unsigned int ooHv = ooLEv + 1;
const unsigned int ooSWv = ooHv + 1;
const unsigned int ooLWv = ooSWv + 1;
const unsigned int ooSWin = ooLWv + 1;
const unsigned int ooLWin = ooSWin + 1;
const unsigned int oomasserror = ooLWin + 1;
const unsigned int ootimestep = oomasserror + 1;
const unsigned int ootot = ootimestep + 1; /*TOTAL NUMBER*/

//****************************************************
//Files
//****************************************************
// first letter "f" ordinary files of input and output
const unsigned int ftsteps = 0;         /*file with time steps*/
const unsigned int fspar = ftsteps + 1; /*soil parameters*/
const unsigned int fmet = fspar + 1;    /*meteo*/
const unsigned int fmetstlist = fmet + 1;
const unsigned int fLRs = fmetstlist + 1; /*lapse rates*/
const unsigned int fhormet = fLRs + 1;    /*horizon of meteo stations*/
const unsigned int fpointlist = fhormet + 1;
const unsigned int fhorpoint = fpointlist + 1; /*horizon of points for which the simulation is run 1D (point_sim==1)*/
const unsigned int fvegpar = fhorpoint + 1; /*vegetation parameter*/
const unsigned int fqin = fvegpar + 1;
const unsigned int fdem = fqin + 1; /*digital elevation model (m)*/
const unsigned int flu = fdem + 1;  /*land use*/
const unsigned int fsoil = flu + 1; /*soil type map*/
const unsigned int fdelay = fsoil + 1;
const unsigned int fsky = fdelay + 1; /*sky view factor*/
const unsigned int fslp = fsky + 1;   /*slope*/
const unsigned int fnet = fslp + 1;   /*channel network*/
const unsigned int fasp = fnet + 1;   /*aspect (0 north, then clockwise)*/
const unsigned int fcurv = fasp + 1;  /*curvature*/
const unsigned int fbed = fcurv + 1;  /*bedrock topography (m)*/
const unsigned int fwt0 = fbed + 1;
const unsigned int fsn0 = fwt0 + 1; /*initial snow depth (mm)*/
const unsigned int fswe0 = fsn0 + 1;
const unsigned int fsnag0 = fswe0 + 1; /*initial snow age (days)*/
const unsigned int fgl0 = fsnag0 + 1;  /*initial glacier depth (mm)*/
const unsigned int fQ = fgl0 + 1;      /*(o.) output discharge file*/

const unsigned int fbas = fQ + 1; /*o. basin variables*/
const unsigned int fbaswriteend = fbas + 1;

const unsigned int fpoint = fbaswriteend + 1; /*o. point variables*/
const unsigned int fpointwriteend = fpoint + 1;

const unsigned int fTz = fpointwriteend + 1; /*o. temperature profiles*/
const unsigned int fTzwriteend = fTz + 1;

const unsigned int fTzav = fTzwriteend + 1;
const unsigned int fTzavwriteend = fTzav + 1;

const unsigned int fpsiz = fTzavwriteend + 1; /*o. psi profiles*/
const unsigned int fpsizwriteend = fpsiz + 1;

const unsigned int fpsiztot = fpsizwriteend + 1; /*o. psi profiles*/
const unsigned int fpsiztotwriteend = fpsiztot + 1;

const unsigned int fliqz = fpsiztotwriteend + 1; /*o. water content profiles*/
const unsigned int fliqzwriteend = fliqz + 1;

const unsigned int fliqzav = fliqzwriteend + 1; /*o. water content profiles*/
const unsigned int fliqzavwriteend = fliqzav + 1;

const unsigned int ficez = fliqzavwriteend + 1; /*o. ice content profiles*/
const unsigned int ficezwriteend = ficez + 1;

const unsigned int ficezav = ficezwriteend + 1; /*o. ice content profiles*/
const unsigned int ficezavwriteend = ficezav + 1;

const unsigned int fsatz = ficezavwriteend + 1;

const unsigned int fsnTz = fsatz + 1; /*o. snow data*/
const unsigned int fsnlz = fsnTz + 1;
const unsigned int fsniz = fsnlz + 1;
const unsigned int fsndz = fsniz + 1;

const unsigned int fsnTzwriteend = fsndz + 1; /*o. snow data*/
const unsigned int fsnlzwriteend = fsnTzwriteend + 1;
const unsigned int fsnizwriteend = fsnlzwriteend + 1;
const unsigned int fsndzwriteend = fsnizwriteend + 1;

const unsigned int fglz = fsndzwriteend + 1; /*o. glacier data*/
const unsigned int fglzwriteend = fglz + 1;  /*o. glacier data*/

const unsigned int fSCA = fglzwriteend + 1; /*file giving the fraction of snow free areas and corresponding properties*/

const unsigned int fTrun = fSCA + 1;
const unsigned int fwrun = fTrun + 1;
const unsigned int fdUrun = fwrun + 1;
const unsigned int fSWErun = fdUrun + 1;
const unsigned int fTmaxrun = fSWErun + 1;
const unsigned int fTminrun = fTmaxrun + 1;
const unsigned int fwmaxrun = fTminrun + 1;
const unsigned int fwminrun = fwmaxrun + 1;

const unsigned int fT = fwminrun + 1; /*o. temperature maps*/
const unsigned int fTsup = fT + 1;    /*o. temperature maps*/
const unsigned int fTav = fTsup + 1;
const unsigned int fTavsup = fTav + 1;
const unsigned int fliq = fTavsup + 1; /*o. water content maps*/
const unsigned int fliqsup = fliq + 1; /*o. water content in the soil at the surface*/
const unsigned int fliqav = fliqsup + 1;
const unsigned int fice = fliqav + 1;  /*o. ice content maps*/
const unsigned int ficesup = fice + 1; /*o. ice content maps*/
const unsigned int ficeav = ficesup + 1;
const unsigned int fhsupland = ficeav + 1; /*o. water over the surface (mm) maps*/
const unsigned int fhsupch = fhsupland + 1; /*o. water over the surface (mm) maps*/

const unsigned int fradnet = fhsupch + 1; /*o. radiation maps*/
const unsigned int fradLWin = fradnet + 1;
const unsigned int fradLW = fradLWin + 1;
const unsigned int fradSW = fradLW + 1;
const unsigned int fradSWin = fradSW + 1;
const unsigned int fradSWinbeam = fradSWin + 1;
const unsigned int fshadow = fradSWinbeam + 1;

const unsigned int fG = fshadow + 1; /*o. surface heat flux maps*/
const unsigned int fH = fG + 1;      /*o. sensible heat flux maps*/
const unsigned int fLE = fH + 1;     /*o. latent heat flux maps*/
const unsigned int fTs = fLE + 1;    /*o. surface temperature maps*/
const unsigned int fprec = fTs + 1;  /*o. precipitation maps*/
const unsigned int fcint = fprec + 1; /*o. precipitation intercepted by canopy maps*/
const unsigned int fpsiliq = fcint + 1; /*o. psi maps*/
const unsigned int fpsitot = fpsiliq + 1;
const unsigned int fsnowdepth = fpsitot + 1;    /*o. snow maps*/
const unsigned int fglacdepth = fsnowdepth + 1; /*o. glacier maps*/
const unsigned int fsnowmelt = fglacdepth + 1;  /*o. snow melted maps*/
const unsigned int fsnowsubl = fsnowmelt + 1;   /*o. snow sublimated maps*/
const unsigned int fglacmelt = fsnowsubl + 1;   /*o. glacier ice melted maps*/
const unsigned int fglacsubl = fglacmelt + 1; /*o. glacier ice sublimated maps*/
const unsigned int fTa = fglacsubl + 1;       /*o. air temperature maps*/
const unsigned int fwspd = fTa + 1;           /*o. wind speed maps*/
const unsigned int fwdir = fwspd + 1;         /*o. wind direction maps*/
const unsigned int frh = fwdir + 1;           /*o. relative humidity maps*/
const unsigned int fswe = frh + 1;            /*o. snow density maps*/
const unsigned int fgwe = fswe + 1;           /*o. glacier ice density maps*/
const unsigned int fsndur = fgwe + 1;         /*o. snow duration maps (hrs)*/
const unsigned int fthawed_up = fsndur + 1;   /*thawed soil depth*/
const unsigned int fthawed_dw = fthawed_up + 1;
const unsigned int fwtable_up = fthawed_dw + 1;
const unsigned int fwtable_dw = fwtable_up + 1; /*water table depth*/
const unsigned int fpnet = fwtable_dw + 1;
const unsigned int fevap = fpnet + 1;

// first letter "p" are special plot files
const unsigned int pG = fevap + 1;
const unsigned int pH = pG + 1;
const unsigned int pLE = pH + 1;
const unsigned int pHg = pLE + 1; /*specific day map plots(p.) sensible heat flux*/
const unsigned int pLEg = pHg + 1; /*p. latent heat flux*/
const unsigned int pHv = pLEg + 1;
const unsigned int pLEv = pHv + 1;
const unsigned int pSWin = pLEv + 1; /*p. incoming shortwave radiation*/
const unsigned int pSWg = pSWin + 1; /*p. outgoing shortwave radiation*/
const unsigned int pSWv = pSWg + 1;
const unsigned int pLWin = pSWv + 1; /*p. incoming longwave radiation*/
const unsigned int pLWg = pLWin + 1; /*p. outgoing longwave radiation*/
const unsigned int pLWv = pLWg + 1;
const unsigned int pTs = pLWv + 1; /*p. surface temperature*/
const unsigned int pTg = pTs + 1;
const unsigned int pTv = pTg + 1;
const unsigned int pTa = pTv + 1;     /*p. air temperature*/
const unsigned int pVspd = pTa + 1;   /*p. wind speed*/
const unsigned int pVdir = pVspd + 1; /*p. wind direction*/
const unsigned int pRH = pVdir + 1;   /*p. relative humidity*/
const unsigned int pD = pRH + 1;      /*p. snow depth*/
const unsigned int pth = pD + 1; /*p. water content of the most superficial layer*/

// first letter "r" are recovery files
const unsigned int rpsi = pth + 1; /*recover file (f.) psi (liquid water pressure)*/
const unsigned int riceg = rpsi + 1;  /*r. soil ice content*/
const unsigned int rTg = riceg + 1;   /*r. soil temperature*/
const unsigned int rDzs = rTg + 1;    /*r. snow layer thicknesses*/
const unsigned int rwls = rDzs + 1;   /*r. snow liquid water contents*/
const unsigned int rwis = rwls + 1;   /*r. snow ice contents*/
const unsigned int rTs = rwis + 1;    /*r. snow temperatures*/
const unsigned int rDzi = rTs + 1;    /*r. glacier layer thicknesses*/
const unsigned int rwli = rDzi + 1;   /*r. glacier liquid water contents*/
const unsigned int rwii = rwli + 1;   /*r. glacier ice contents*/
const unsigned int rTi = rwii + 1;    /*r. glacier temperatures*/
const unsigned int rns = rTi + 1;     /*r. number of snow layers*/
const unsigned int rni = rns + 1;     /*r. number of glacier layers*/
const unsigned int rsnag = rni + 1;   /*r. snow age*/
const unsigned int rwcrn = rsnag + 1; /*r. liquid water stored in canopy*/
const unsigned int rwcsn = rwcrn + 1; /*r. snow stored in canopy*/
const unsigned int rTv = rwcsn + 1;   /*r. canopy temperature*/
const unsigned int rpsich = rTv + 1; /*r. psi (liquid water pressure) in channels*/
const unsigned int ricegch = rpsich + 1;
const unsigned int rTgch = ricegch + 1;
const unsigned int rTrun = rTgch + 1;
const unsigned int rwrun = rTrun + 1;
const unsigned int rdUrun = rwrun + 1;
const unsigned int rSWErun = rdUrun + 1;
const unsigned int rTmaxrun = rSWErun + 1;
const unsigned int rTminrun = rTmaxrun + 1;
const unsigned int rwmaxrun = rTminrun + 1;
const unsigned int rwminrun = rwmaxrun + 1;
const unsigned int rtime = rwminrun + 1;
const unsigned int rsux = rtime + 1;

const unsigned int nfiles = rsux + 1; /*number of files*/

//****************************************************
//Points
//****************************************************
const unsigned int ptID = 1;
const unsigned int ptX = ptID + 1;
const unsigned int ptY = ptX + 1;
const unsigned int ptZ = ptY + 1;

const unsigned int ptLC = ptZ + 1;
const unsigned int ptSY = ptLC + 1;
const unsigned int ptS = ptSY + 1;
const unsigned int ptA = ptS + 1;
const unsigned int ptSKY = ptA + 1;
const unsigned int ptCNS = ptSKY + 1;
const unsigned int ptCWE = ptCNS + 1;
const unsigned int ptCNwSe = ptCWE + 1;
const unsigned int ptCNeSw = ptCNwSe + 1;
const unsigned int ptDrDEPTH = ptCNeSw + 1;
const unsigned int ptHOR = ptDrDEPTH + 1;
const unsigned int ptMAXSWE = ptHOR + 1;
const unsigned int ptLAT = ptMAXSWE + 1;
const unsigned int ptLON = ptLAT + 1;
const unsigned int ptBED = ptLON + 1;
const unsigned int ptTOT = ptBED;

#endif
