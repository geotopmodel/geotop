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
    constexpr double KNe = 0.0;  // Euler method parameter for heat equation (0.5 = Crank Nicholson, 0 = Backward Euler)
    constexpr double LSAIthres = 0.1;  // Minimum LSAI
    constexpr double z_evap = 100.0;  // soil depth responsable for soil evaporation [mm]
    constexpr double z_transp = 10000.0;  // soil depth responsable for canopy transpiration [mm]
    constexpr double min_tau_cloud = 0.1;
    constexpr double RelativeErrorRichards = 1.E-10;
    constexpr int max_cols_time_steps_file = 100;
    constexpr double PsiMin = -1.E+10;
    constexpr double thmin = 0.1;  // Newton method parameter
    constexpr double thmax = 0.5;  // Newton method parameter
    constexpr double Tol_h_mount = 2.0;
    constexpr double Tol_h_flat = 12.0;
    constexpr double max_slope = 89.999;

    // STANDARD LAPSE RATES: For the sign, remember that the Lapse Rate gives how a variable decrease with height
    constexpr double LapseRateTair = 6.5;  /** Lapse rate for Tair [°C/m] */
    constexpr double LapseRateTdew = 2.5;  /** Lapse rate for Tdew [°C/m] */
    constexpr double LapseRatePrec = 0.0;  /** Lapse rate for Precipitation [1/m] */
//****************************************************
// Constants
//****************************************************
    constexpr double omega = 0.261799388; /* Earth Rotation Velocity [rad/h] */
    constexpr double Isc = 1367;          /* Solar Constant  [W/m2] */
    constexpr double Pa0 = 1013.25;       /* Mean atmospheric at sea level [mbar] */
    constexpr double rho_w = 1000;        /* density of water [kg/m3] */
    constexpr double rho_i = 917;         /* density of ice [kg/m3] */
    constexpr double Lf = 333700.00;      /* latent heat of fusion [J/kg] */
    constexpr double Ls = 2834700.00;     /* latent heat of sublimation [J/kg] */
    constexpr double g = 9.81;      /* gravity acceleration [m/s2] */
    constexpr double Pi = 3.14159265358979;
    constexpr double tk = 273.15;    /* =0 Deg in Kelvin*/
    constexpr double k_liq = 0.567;  /* thermal conductivity of water [W m^-1 K^-1]*/
    constexpr double k_ice = 2.290;  /* thermal conductivity of water [W m^-1 K^-1]*/
    constexpr double k_air = 0.023;  /* thermal conductivity of air   [W m^-1 K^-1]*/
    constexpr double c_liq = 4188.0; /* heat capacity of water    [J/(kg*K)]*/
    constexpr double c_ice = 2117.0; /* heat capacity of ice    [J/(kg*K)]*/
    constexpr double c_can = 2700.0; /* heat capacity of canopy [J/(kg*K)]*/
    constexpr double Tfreezing = 0.0; /** freezing temperature [*C] */
    constexpr double ka = 0.41;         /* Von Karman constant */
    constexpr double mu_l = 0.001787; /* Dynamic viscosity of water at 0 degrees Celsius*/
    constexpr double wsn_vis = 0.8; /* snow on canopy: scattering parameters */
    constexpr double wsn_nir = 0.4;
    constexpr double Bsnd_vis = 0.5;
    constexpr double Bsnd_nir = 0.5;
    constexpr double Bsnb_vis = 0.5;
    constexpr double Bsnb_nir = 0.5;
    constexpr double Cd = 0.61; /* discharge coefficient of the weir model */
    constexpr double Rwv = 461.495;  // Specific gas constant for water vapor, 461.495 J/(kg·K)
    constexpr double D00 = 21.7;  // molecular diffusivity of water vapor, 21.7 mm2/s
    constexpr double secinday = 86400.0;  // seconds in one day
} // end namespace GTConst

//****************************************************
//Meteo data
//****************************************************
constexpr unsigned int iDate12 = 0;             /*Date12 : DDMMYYYYhhmm*/
constexpr unsigned int iJDfrom0 = iDate12 + 1;  /*Julian Day from year 0*/
constexpr unsigned int iPrecInt = iJDfrom0 + 1; /*Precipitation*/
constexpr unsigned int iPrec = iPrecInt + 1;
constexpr unsigned int iWs = iPrec + 1;  /*Total wind speed*/
constexpr unsigned int iWdir = iWs + 1;  /*Wind direction*/
constexpr unsigned int iWsx = iWdir + 1; /*Wind speed from west, to east*/
constexpr unsigned int iWsy = iWsx + 1;  /*Wind speed from south, to north*/
constexpr unsigned int iRh = iWsy + 1;   /*Relative humidity*/
constexpr unsigned int iT = iRh + 1;     /*Air temperature*/
constexpr unsigned int iTdew = iT + 1;   /*Air dew temperature*/
constexpr unsigned int iSW = iTdew + 1;  /*Global shortwave radiation*/
constexpr unsigned int iSWb = iSW + 1;   /*Direct SW*/
constexpr unsigned int iSWd = iSWb + 1;  /*Diffuse SW*/
constexpr unsigned int itauC = iSWd + 1; /*Cloud transmissivity in SWin*/
constexpr unsigned int iC = itauC + 1;   /*Cloudiness factor*/
constexpr unsigned int iLWi = iC + 1;    /*Incoming longwave*/
constexpr unsigned int iSWn = iLWi + 1;  /*Net shortwave*/
constexpr unsigned int iTs = iSWn + 1;   /*Surface Temperature*/
constexpr unsigned int iTbottom = iTs + 1;   /*Bottom Temperature (NOT present in GEOtop 2.1)*/
constexpr unsigned int nmet = iTbottom + 1;

//****************************************************
//soil data
//****************************************************
constexpr unsigned int ilsDate12 = 0;
constexpr unsigned int ilsTa = ilsDate12 + 1;
constexpr unsigned int ilsTdew = ilsTa + 1;
constexpr unsigned int ilsPrec = ilsTdew + 1;
constexpr unsigned int nlstot = ilsPrec + 1;

//****************************************************
//soil data
//****************************************************
constexpr unsigned int jdz = 1;         /*layer thickness [mm]*/
constexpr unsigned int jpsi = jdz + 1;  /*initial psi [mm]*/
constexpr unsigned int jT = jpsi + 1;   /*initial temperature [C]*/
constexpr unsigned int jKn = jT + 1;    /*normal hydr. conductivity [mm/s]*/
constexpr unsigned int jKl = jKn + 1;   /*lateral hydr. conductivity [mm/s]*/
constexpr unsigned int jres = jKl + 1;  /*residual wat.cont.*/
constexpr unsigned int jwp = jres + 1;  /*wilting point water cont.*/
constexpr unsigned int jfc = jwp + 1;   /*field capacity water cont. */
constexpr unsigned int jsat = jfc + 1;  /*porosity*/
constexpr unsigned int ja = jsat + 1;   /*alpha[mm^-1]*/
constexpr unsigned int jns = ja + 1;    /*n*/
constexpr unsigned int jv = jns + 1;    /*v*/
constexpr unsigned int jkt = jv + 1;    /*thermal conductivity*/
constexpr unsigned int jct = jkt + 1;   /*thermal capacity*/
constexpr unsigned int jss = jct + 1;   /*soil specific storativity*/
constexpr unsigned int nsoilprop = jss; /*number of soil properties considered*/

//****************************************************
//land use data
//****************************************************
constexpr unsigned int jz0 = 1;     /*roughness length for soil*/
constexpr unsigned int jz0thressoil = jz0 + 1;   /*threshold on snow depth to change roughness length to snow covered values in soil area*/
constexpr unsigned int jHveg = jz0thressoil + 1; /*vegetation height*/
constexpr unsigned int jz0thresveg = jHveg + 1; /*threshold on snow depth to change roughness length to snow covered values in vegetated area*/
constexpr unsigned int jz0thresveg2 = jz0thresveg + 1;
constexpr unsigned int jLSAI = jz0thresveg2 + 1; /*LSAI*/
constexpr unsigned int jcf = jLSAI + 1;          /*Canopy fraction*/
constexpr unsigned int jdecay0 = jcf + 1;
constexpr unsigned int jexpveg = jdecay0 + 1;
constexpr unsigned int jroot = jexpveg + 1; /*root depth [mm]*/
constexpr unsigned int jrs = jroot + 1;     /*canopy transpiration coefficient*/
constexpr unsigned int jvR_vis = jrs + 1;   /*vegetation in the visible spectrum*/
constexpr unsigned int jvR_nir = jvR_vis + 1; /*vegetation in the near infrared spectrum*/
constexpr unsigned int jvT_vis = jvR_nir + 1; /*vegetation in the visible spectrum*/
constexpr unsigned int jvT_nir = jvT_vis + 1; /*vegetation in the near infrared spectrum*/
constexpr unsigned int jvCh = jvT_nir + 1; /*departure of leaf angles from a random distribution (1 horizontal, 0 random, -1 vertical)*/
constexpr unsigned int jcd = jvCh + 1; /*surface density of canopy [kg/(m2*LSAI)]*/
constexpr unsigned int ja_vis_dry = jcd + 1; /*ground albedo in the visible spectrum dry*/
constexpr unsigned int ja_nir_dry = ja_vis_dry + 1; /*ground albedo in the near infrared spectrum dry*/
constexpr unsigned int ja_vis_sat = ja_nir_dry + 1; /*ground albedo in the visible spectrum saturated*/
constexpr unsigned int ja_nir_sat = ja_vis_sat + 1; /*ground albedo in the near infrared spectrum saturated*/
constexpr unsigned int jemg = ja_nir_sat + 1; /*soil emissivity*/
constexpr unsigned int jcm = jemg + 1; /*gauckler strickler (1/manning) coefficient*/
constexpr unsigned int jN = jcm + 1; /*number of roughness elements (vegetation) per unit surface for blowing snow*/
constexpr unsigned int jdv = jN + 1; /*diameter of roughness elements (vegetation) [mm]*/
constexpr unsigned int nlandprop = jdv; /*number of land use properties*/

//****************************************************
//vegetation files (the numbers must respect the same order as the block above)
//****************************************************
constexpr unsigned int jdHveg = 1; /*vegetation height*/
constexpr unsigned int jdz0thresveg = jdHveg + 1; /*threshold on snow depth to change roughness length to snow covered values in vegetated area*/
constexpr unsigned int jdz0thresveg2 = jdz0thresveg + 1;
constexpr unsigned int jdLSAI = jdz0thresveg2 + 1; /*LSAI*/
constexpr unsigned int jdcf = jdLSAI + 1;          /*Canopy fraction*/
constexpr unsigned int jddecay0 = jdcf + 1;
constexpr unsigned int jdexpveg = jddecay0 + 1;
constexpr unsigned int jdroot = jdexpveg + 1;
constexpr unsigned int jdrs = jdroot + 1;
constexpr unsigned int jdvegprop = jdrs;

//****************************************************
//point output
//****************************************************
constexpr unsigned int odate12 = 0;
constexpr unsigned int oJDfrom0 = odate12 + 1;
constexpr unsigned int odaysfromstart = oJDfrom0 + 1;
constexpr unsigned int operiod = odaysfromstart + 1;
constexpr unsigned int orun = operiod + 1;
constexpr unsigned int opoint = orun + 1;
constexpr unsigned int osnowover = opoint + 1;    /*prec_snow_atm;*/
constexpr unsigned int orainover = osnowover + 1; /*prec_rain_atm;*/
constexpr unsigned int oprecsnow = orainover + 1; /*prec_snow;*/
constexpr unsigned int oprecrain = oprecsnow + 1; /*(prec_rain_on_soil+prec_rain_on_snow);*/
constexpr unsigned int orainonsnow = oprecrain + 1; /*prec_rain_on_snow;*/
constexpr unsigned int oV = orainonsnow + 1;        /*Vpoint/(double)n;*/
constexpr unsigned int oVdir = oV + 1;   /*(met->Vdir->co[r][c])/(double)n;*/
constexpr unsigned int oRH = oVdir + 1;  /*RHpoint/(double)n;*/
constexpr unsigned int oPa = oRH + 1;    /*Ppoint/(double)n;*/
constexpr unsigned int oTa = oPa + 1;    /*Tpoint/(double)n;*/
constexpr unsigned int oTdew = oTa + 1;  /*Tdew/(double)n;*/
constexpr unsigned int oTg = oTdew + 1;  /*Tg/(double)n; //Ts[C]*/
constexpr unsigned int oTv = oTg + 1;    /*Tv/(double)n;*/
constexpr unsigned int oTs = oTv + 1;    /*Ts/(double)n;*/
constexpr unsigned int oEB = oTs + 1;    /*surfEB/(double)n;*/
constexpr unsigned int oG = oEB + 1;     /*G/(double)n;*/
constexpr unsigned int oSWin = oG + 1;   /*SWin/(double)n;*/
constexpr unsigned int oSWb = oSWin + 1; /*(SWbeam/(double)n);*/
constexpr unsigned int oSWd = oSWb + 1;  /*(SWdiff/(double)n);*/
constexpr unsigned int oLWin = oSWd + 1; /*LWin/(double)n;*/
constexpr unsigned int ominLWin = oLWin + 1; /*(epsa_min*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;*/
constexpr unsigned int omaxLWin = ominLWin + 1; /*(epsa_max*5.67E-8*pow(Tpoint+tk,4.0))/(double)n;*/
constexpr unsigned int oSW = omaxLWin + 1;   /*SW/(double)n;*/
constexpr unsigned int oLW = oSW + 1;        /*LW/(double)n;*/
constexpr unsigned int oH = oLW + 1;         /*H/(double)n; //H[W/m^2]*/
constexpr unsigned int oLE = oH + 1;         /*LE/(double)n; //ET[W/m^2]*/
constexpr unsigned int ofc = oLE + 1;        /*fc/(double)n;*/
constexpr unsigned int oLSAI = ofc + 1;      /*LSAI/(double)n;*/
constexpr unsigned int oz0v = oLSAI + 1;     /*z0/(double)n;*/
constexpr unsigned int od0v = oz0v + 1;      /*d0/(double)n;*/
constexpr unsigned int oEcan = od0v + 1;     /*(SWv+LWv-Hv-LEv)/(double)n;*/
constexpr unsigned int oSWv = oEcan + 1;     /*SWv/(double)n;*/
constexpr unsigned int oLWv = oSWv + 1;      /*LWv/(double)n;*/
constexpr unsigned int oHv = oLWv + 1;       /*Hv/(double)n;*/
constexpr unsigned int oLEv = oHv + 1;       /*LEv/(double)n;*/
constexpr unsigned int oHg0 = oLEv + 1;      /*Hg0/(double)n;*/
constexpr unsigned int oLEg0 = oHg0 + 1;     /*Levap(Tg)*Eg0/(double)n;*/
constexpr unsigned int oHg1 = oLEg0 + 1;     /*Hg1/(double)n;*/
constexpr unsigned int oLEg1 = oHg1 + 1;     /*Levap(Tg)*Eg1/(double)n;*/
constexpr unsigned int oevapsur = oLEg1 + 1; /*Er_soil*par->Dt; //Eg[mm]*/
constexpr unsigned int otrasp = oevapsur + 1; /*Evt*(1000.0/rho_w)*par->Dt; //Etc[mm]*/
constexpr unsigned int owcan_rain = otrasp + 1; /*wat->wcan_rain->co[r][c]/(double)n;*/
constexpr unsigned int owcan_snow = owcan_rain + 1;                        /*wat->wcan_snow->co[r][c]/(double)n;*/
constexpr unsigned int oQv = owcan_snow + 1; /*(Qv)/(double)n;*/
constexpr unsigned int oQg = oQv + 1;        /*(Qg)/(double)n;*/
constexpr unsigned int oQa = oQg + 1;        /*(Qa)/(double)n;*/
constexpr unsigned int oQs = oQa + 1;        /*(Qs)/(double)n;*/
constexpr unsigned int oLobuk = oQs + 1;     /*turbulence->co[2]/(double)n;*/
constexpr unsigned int oLobukcan = oLobuk + 1; /*(Locc)/(double)n;*/
constexpr unsigned int outop = oLobukcan + 1;  /*(u_top)/(double)n;*/
constexpr unsigned int odecay = outop + 1;     /*(decay)/(double)n;*/
constexpr unsigned int oSWup = odecay + 1;     /*SWupabove_v/(double)n;*/
constexpr unsigned int oLWup = oSWup + 1;      /*LWup_above_v/(double)n;*/
constexpr unsigned int oHup = oLWup + 1;       /*(H+fc*Hv)/(double)n;*/
constexpr unsigned int oLEup = oHup + 1;       /*(LE+fc*LEv)/(double)n;*/
constexpr unsigned int osnowdepth = oLEup + 1;
constexpr unsigned int oSWE = osnowdepth + 1;
constexpr unsigned int osnowdens = oSWE + 1;
constexpr unsigned int osnowT = osnowdens + 1;
constexpr unsigned int omrsnow = osnowT + 1;  /*Mr_snow*par->Dt;  //[mm]*/
constexpr unsigned int osrsnow = omrsnow + 1; /*Sr_snow*par->Dt;  //[mm]*/
constexpr unsigned int oblowingsnowtrans = osrsnow + 1;
constexpr unsigned int oblowingsnowsubl = oblowingsnowtrans + 1;
constexpr unsigned int oglacdepth = oblowingsnowsubl + 1;
constexpr unsigned int oGWE = oglacdepth + 1;
constexpr unsigned int oglacdens = oGWE + 1;
constexpr unsigned int oglacT = oglacdens + 1;
constexpr unsigned int omrglac = oglacT + 1;      /*Mr_glac*par->Dt;  //[mm]*/
constexpr unsigned int osrglac = omrglac + 1;     /*Sr_glac*par->Dt;  //[mm]*/
constexpr unsigned int othawedup = osrglac + 1;   /*thawed soil depth [mm]*/
constexpr unsigned int othaweddw = othawedup + 1; /*thawed soil depth [mm]*/
constexpr unsigned int owtableup = othaweddw + 1; /*water table depth [mm]*/
constexpr unsigned int owtabledw = owtableup + 1; /*water table depth [mm]*/
constexpr unsigned int otot = owtabledw + 1;      /*TOTAL NUMBER*/

//****************************************************
//BASIN OUTPUT
//****************************************************
constexpr unsigned int oodate12 = 0;
constexpr unsigned int ooJDfrom0 = oodate12 + 1;
constexpr unsigned int oodaysfromstart = ooJDfrom0 + 1;
constexpr unsigned int ooperiod = oodaysfromstart + 1;
constexpr unsigned int oorun = ooperiod + 1;
constexpr unsigned int ooprecrain = oorun + 1;
constexpr unsigned int ooprecsnow = ooprecrain + 1;
constexpr unsigned int oorainover = ooprecsnow + 1;
constexpr unsigned int oosnowover = oorainover + 1;
constexpr unsigned int oopnet = oosnowover + 1;
constexpr unsigned int ooTa = oopnet + 1;
constexpr unsigned int ooTg = ooTa + 1;
constexpr unsigned int ooTv = ooTg + 1;
constexpr unsigned int ooevapsur = ooTv + 1;
constexpr unsigned int ootrasp = ooevapsur + 1;
constexpr unsigned int ooLE = ootrasp + 1;
constexpr unsigned int ooH = ooLE + 1;
constexpr unsigned int ooSW = ooH + 1;
constexpr unsigned int ooLW = ooSW + 1;
constexpr unsigned int ooLEv = ooLW + 1;
constexpr unsigned int ooHv = ooLEv + 1;
constexpr unsigned int ooSWv = ooHv + 1;
constexpr unsigned int ooLWv = ooSWv + 1;
constexpr unsigned int ooSWin = ooLWv + 1;
constexpr unsigned int ooLWin = ooSWin + 1;
constexpr unsigned int oomasserror = ooLWin + 1;
constexpr unsigned int ootimestep = oomasserror + 1;
constexpr unsigned int ootot = ootimestep + 1; /*TOTAL NUMBER*/

//****************************************************
//Files
//****************************************************
// first letter "f" ordinary files of input and output
constexpr unsigned int ftsteps = 0;         /*file with time steps*/
constexpr unsigned int fspar = ftsteps + 1; /*soil parameters*/
constexpr unsigned int fmet = fspar + 1;    /*meteo*/
constexpr unsigned int fmetstlist = fmet + 1;
constexpr unsigned int fLRs = fmetstlist + 1; /*lapse rates*/
constexpr unsigned int fhormet = fLRs + 1;    /*horizon of meteo stations*/
constexpr unsigned int fpointlist = fhormet + 1;
constexpr unsigned int fhorpoint = fpointlist + 1; /*horizon of points for which the simulation is run 1D (point_sim==1)*/
constexpr unsigned int fvegpar = fhorpoint + 1; /*vegetation parameter*/
constexpr unsigned int fqin = fvegpar + 1;
constexpr unsigned int fdem = fqin + 1; /*digital elevation model (m)*/
constexpr unsigned int flu = fdem + 1;  /*land use*/
constexpr unsigned int fsoil = flu + 1; /*soil type map*/
constexpr unsigned int fdelay = fsoil + 1;
constexpr unsigned int fsky = fdelay + 1; /*sky view factor*/
constexpr unsigned int fslp = fsky + 1;   /*slope*/
constexpr unsigned int fnet = fslp + 1;   /*channel network*/
constexpr unsigned int fasp = fnet + 1;   /*aspect (0 north, then clockwise)*/
constexpr unsigned int fcurv = fasp + 1;  /*curvature*/
constexpr unsigned int fbed = fcurv + 1;  /*bedrock topography (m)*/
constexpr unsigned int fwt0 = fbed + 1;
constexpr unsigned int fsn0 = fwt0 + 1; /*initial snow depth (mm)*/
constexpr unsigned int fswe0 = fsn0 + 1;
constexpr unsigned int fsnag0 = fswe0 + 1; /*initial snow age (days)*/
constexpr unsigned int fgl0 = fsnag0 + 1;  /*initial glacier depth (mm)*/
constexpr unsigned int fQ = fgl0 + 1;      /*(o.) output discharge file*/

constexpr unsigned int fbas = fQ + 1; /*o. basin variables*/
constexpr unsigned int fbaswriteend = fbas + 1;

constexpr unsigned int fpoint = fbaswriteend + 1; /*o. point variables*/
constexpr unsigned int fpointwriteend = fpoint + 1;

constexpr unsigned int fTz = fpointwriteend + 1; /*o. temperature profiles*/
constexpr unsigned int fTzwriteend = fTz + 1;

constexpr unsigned int fTzav = fTzwriteend + 1;
constexpr unsigned int fTzavwriteend = fTzav + 1;

constexpr unsigned int fpsiz = fTzavwriteend + 1; /*o. psi profiles*/
constexpr unsigned int fpsizwriteend = fpsiz + 1;

constexpr unsigned int fpsiztot = fpsizwriteend + 1; /*o. psi profiles*/
constexpr unsigned int fpsiztotwriteend = fpsiztot + 1;

constexpr unsigned int fliqz = fpsiztotwriteend + 1; /*o. water content profiles*/
constexpr unsigned int fliqzwriteend = fliqz + 1;

constexpr unsigned int fliqzav = fliqzwriteend + 1; /*o. water content profiles*/
constexpr unsigned int fliqzavwriteend = fliqzav + 1;

constexpr unsigned int ficez = fliqzavwriteend + 1; /*o. ice content profiles*/
constexpr unsigned int ficezwriteend = ficez + 1;

constexpr unsigned int ficezav = ficezwriteend + 1; /*o. ice content profiles*/
constexpr unsigned int ficezavwriteend = ficezav + 1;

constexpr unsigned int fsatz = ficezavwriteend + 1;

constexpr unsigned int fsnTz = fsatz + 1; /*o. snow data*/
constexpr unsigned int fsnlz = fsnTz + 1;
constexpr unsigned int fsniz = fsnlz + 1;
constexpr unsigned int fsndz = fsniz + 1;

constexpr unsigned int fsnTzwriteend = fsndz + 1; /*o. snow data*/
constexpr unsigned int fsnlzwriteend = fsnTzwriteend + 1;
constexpr unsigned int fsnizwriteend = fsnlzwriteend + 1;
constexpr unsigned int fsndzwriteend = fsnizwriteend + 1;

constexpr unsigned int fglz = fsndzwriteend + 1; /*o. glacier data*/
constexpr unsigned int fglzwriteend = fglz + 1;  /*o. glacier data*/

constexpr unsigned int fSCA = fglzwriteend + 1; /*file giving the fraction of snow free areas and corresponding properties*/

constexpr unsigned int fTrun = fSCA + 1;
constexpr unsigned int fwrun = fTrun + 1;
constexpr unsigned int fdUrun = fwrun + 1;
constexpr unsigned int fSWErun = fdUrun + 1;
constexpr unsigned int fTmaxrun = fSWErun + 1;
constexpr unsigned int fTminrun = fTmaxrun + 1;
constexpr unsigned int fwmaxrun = fTminrun + 1;
constexpr unsigned int fwminrun = fwmaxrun + 1;

constexpr unsigned int fT = fwminrun + 1; /*o. temperature maps*/
constexpr unsigned int fTsup = fT + 1;    /*o. temperature maps*/
constexpr unsigned int fTav = fTsup + 1;
constexpr unsigned int fTavsup = fTav + 1;
constexpr unsigned int fliq = fTavsup + 1; /*o. water content maps*/
constexpr unsigned int fliqsup = fliq + 1; /*o. water content in the soil at the surface*/
constexpr unsigned int fliqav = fliqsup + 1;
constexpr unsigned int fice = fliqav + 1;  /*o. ice content maps*/
constexpr unsigned int ficesup = fice + 1; /*o. ice content maps*/
constexpr unsigned int ficeav = ficesup + 1;
constexpr unsigned int fhsupland = ficeav + 1; /*o. water over the surface (mm) maps*/
constexpr unsigned int fhsupch = fhsupland + 1; /*o. water over the surface (mm) maps*/

constexpr unsigned int fradnet = fhsupch + 1; /*o. radiation maps*/
constexpr unsigned int fradLWin = fradnet + 1;
constexpr unsigned int fradLW = fradLWin + 1;
constexpr unsigned int fradSW = fradLW + 1;
constexpr unsigned int fradSWin = fradSW + 1;
constexpr unsigned int fradSWinbeam = fradSWin + 1;
constexpr unsigned int fshadow = fradSWinbeam + 1;

constexpr unsigned int fG = fshadow + 1; /*o. surface heat flux maps*/
constexpr unsigned int fH = fG + 1;      /*o. sensible heat flux maps*/
constexpr unsigned int fLE = fH + 1;     /*o. latent heat flux maps*/
constexpr unsigned int fTs = fLE + 1;    /*o. surface temperature maps*/
constexpr unsigned int fprec = fTs + 1;  /*o. precipitation maps*/
constexpr unsigned int fcint = fprec + 1; /*o. precipitation intercepted by canopy maps*/
constexpr unsigned int fpsiliq = fcint + 1; /*o. psi maps*/
constexpr unsigned int fpsitot = fpsiliq + 1;
constexpr unsigned int fsnowdepth = fpsitot + 1;    /*o. snow maps*/
constexpr unsigned int fglacdepth = fsnowdepth + 1; /*o. glacier maps*/
constexpr unsigned int fsnowmelt = fglacdepth + 1;  /*o. snow melted maps*/
constexpr unsigned int fsnowsubl = fsnowmelt + 1;   /*o. snow sublimated maps*/
constexpr unsigned int fglacmelt = fsnowsubl + 1;   /*o. glacier ice melted maps*/
constexpr unsigned int fglacsubl = fglacmelt + 1; /*o. glacier ice sublimated maps*/
constexpr unsigned int fTa = fglacsubl + 1;       /*o. air temperature maps*/
constexpr unsigned int fwspd = fTa + 1;           /*o. wind speed maps*/
constexpr unsigned int fwdir = fwspd + 1;         /*o. wind direction maps*/
constexpr unsigned int frh = fwdir + 1;           /*o. relative humidity maps*/
constexpr unsigned int fswe = frh + 1;            /*o. snow density maps*/
constexpr unsigned int fgwe = fswe + 1;           /*o. glacier ice density maps*/
constexpr unsigned int fsndur = fgwe + 1;         /*o. snow duration maps (hrs)*/
constexpr unsigned int fthawed_up = fsndur + 1;   /*thawed soil depth*/
constexpr unsigned int fthawed_dw = fthawed_up + 1;
constexpr unsigned int fwtable_up = fthawed_dw + 1;
constexpr unsigned int fwtable_dw = fwtable_up + 1; /*water table depth*/
constexpr unsigned int fpnet = fwtable_dw + 1;
constexpr unsigned int fevap = fpnet + 1;

// first letter "p" are special plot files
constexpr unsigned int pG = fevap + 1;
constexpr unsigned int pH = pG + 1;
constexpr unsigned int pLE = pH + 1;
constexpr unsigned int pHg = pLE + 1; /*specific day map plots(p.) sensible heat flux*/
constexpr unsigned int pLEg = pHg + 1; /*p. latent heat flux*/
constexpr unsigned int pHv = pLEg + 1;
constexpr unsigned int pLEv = pHv + 1;
constexpr unsigned int pSWin = pLEv + 1; /*p. incoming shortwave radiation*/
constexpr unsigned int pSWg = pSWin + 1; /*p. outgoing shortwave radiation*/
constexpr unsigned int pSWv = pSWg + 1;
constexpr unsigned int pLWin = pSWv + 1; /*p. incoming longwave radiation*/
constexpr unsigned int pLWg = pLWin + 1; /*p. outgoing longwave radiation*/
constexpr unsigned int pLWv = pLWg + 1;
constexpr unsigned int pTs = pLWv + 1; /*p. surface temperature*/
constexpr unsigned int pTg = pTs + 1;
constexpr unsigned int pTv = pTg + 1;
constexpr unsigned int pTa = pTv + 1;     /*p. air temperature*/
constexpr unsigned int pVspd = pTa + 1;   /*p. wind speed*/
constexpr unsigned int pVdir = pVspd + 1; /*p. wind direction*/
constexpr unsigned int pRH = pVdir + 1;   /*p. relative humidity*/
constexpr unsigned int pD = pRH + 1;      /*p. snow depth*/
constexpr unsigned int pth = pD + 1; /*p. water content of the most superficial layer*/

// first letter "r" are recovery files
constexpr unsigned int rpsi = pth + 1; /*recover file (f.) psi (liquid water pressure)*/
constexpr unsigned int riceg = rpsi + 1;  /*r. soil ice content*/
constexpr unsigned int rTg = riceg + 1;   /*r. soil temperature*/
constexpr unsigned int rDzs = rTg + 1;    /*r. snow layer thicknesses*/
constexpr unsigned int rwls = rDzs + 1;   /*r. snow liquid water contents*/
constexpr unsigned int rwis = rwls + 1;   /*r. snow ice contents*/
constexpr unsigned int rTs = rwis + 1;    /*r. snow temperatures*/
constexpr unsigned int rDzi = rTs + 1;    /*r. glacier layer thicknesses*/
constexpr unsigned int rwli = rDzi + 1;   /*r. glacier liquid water contents*/
constexpr unsigned int rwii = rwli + 1;   /*r. glacier ice contents*/
constexpr unsigned int rTi = rwii + 1;    /*r. glacier temperatures*/
constexpr unsigned int rns = rTi + 1;     /*r. number of snow layers*/
constexpr unsigned int rni = rns + 1;     /*r. number of glacier layers*/
constexpr unsigned int rsnag = rni + 1;   /*r. snow age*/
constexpr unsigned int rwcrn = rsnag + 1; /*r. liquid water stored in canopy*/
constexpr unsigned int rwcsn = rwcrn + 1; /*r. snow stored in canopy*/
constexpr unsigned int rTv = rwcsn + 1;   /*r. canopy temperature*/
constexpr unsigned int rpsich = rTv + 1; /*r. psi (liquid water pressure) in channels*/
constexpr unsigned int ricegch = rpsich + 1;
constexpr unsigned int rTgch = ricegch + 1;
constexpr unsigned int rTrun = rTgch + 1;
constexpr unsigned int rwrun = rTrun + 1;
constexpr unsigned int rdUrun = rwrun + 1;
constexpr unsigned int rSWErun = rdUrun + 1;
constexpr unsigned int rTmaxrun = rSWErun + 1;
constexpr unsigned int rTminrun = rTmaxrun + 1;
constexpr unsigned int rwmaxrun = rTminrun + 1;
constexpr unsigned int rwminrun = rwmaxrun + 1;
constexpr unsigned int rtime = rwminrun + 1;
constexpr unsigned int rsux = rtime + 1;

constexpr unsigned int nfiles = rsux + 1; /*number of files*/

//****************************************************
//Points
//****************************************************
constexpr unsigned int ptID = 1;
constexpr unsigned int ptX = ptID + 1;
constexpr unsigned int ptY = ptX + 1;
constexpr unsigned int ptZ = ptY + 1;

constexpr unsigned int ptLC = ptZ + 1;
constexpr unsigned int ptSY = ptLC + 1;
constexpr unsigned int ptS = ptSY + 1;
constexpr unsigned int ptA = ptS + 1;
constexpr unsigned int ptSKY = ptA + 1;
constexpr unsigned int ptCNS = ptSKY + 1;
constexpr unsigned int ptCWE = ptCNS + 1;
constexpr unsigned int ptCNwSe = ptCWE + 1;
constexpr unsigned int ptCNeSw = ptCNwSe + 1;
constexpr unsigned int ptDrDEPTH = ptCNeSw + 1;
constexpr unsigned int ptHOR = ptDrDEPTH + 1;
constexpr unsigned int ptMAXSWE = ptHOR + 1;
constexpr unsigned int ptLAT = ptMAXSWE + 1;
constexpr unsigned int ptLON = ptLAT + 1;
constexpr unsigned int ptBED = ptLON + 1;
constexpr unsigned int ptTOT = ptBED;

#endif
