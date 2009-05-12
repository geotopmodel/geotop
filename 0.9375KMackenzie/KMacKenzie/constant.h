
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


/* Name of the program */

#define PROGRAM_NAME  "___geotop"

//Constants
#define omega 0.261799388			/* velocita' di rotazione terrestre [rad/hr] */
#define omega_anno 0.0007167		/* velocita' di rivoluzione terrestre [rad/hr] */
#define Isc 1367					/* Costante solare [W/mq] */
#define Pa0	1013.25					/* Mean atmospheric at sea level [mbar] */
#define rho_w 1000					/* density of water [kg/mc] */
#define rho_i 917					/* density of ice [kg/mc] */
#define Lf 333700.00				/* heat of fusion [J/kg] */
#define g 9.81						/* gravity acceleration [m/s2] */
#define Pi 3.14159265358979			/* greek P*/
#define k_w 0.0000009				/* water heat diffusivity [m^2/s]*/
#define tk 273.15					/* =0 Deg in Kelvin*/
#define k_liq 0.600					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_ice 2.290					/* thermal conductivity of water [W m^-1 K^-1]*/
#define k_air 0.023					/* thermal conductivity of air   [W m^-1 K^-1]*/
#define c_liq 4188.0				/* heat capacity of water		[J/(kg*K)]*/
#define c_ice 2117.0				/* heat capacity of ice		[J/(kg*K)]*/
#define c_can 2700.0				/* heat capacity of canopy [J/(kg*K)]*/
#define CA 0.34						/* tunable parameter in order to have surface temperature instead layer center temperature*/
#define KNe 0.50					/* Krank-Nicholson parameter for egy balance*/
#define KNw 1.00					/* Krank-Nicholson parameter for water balance*/
#define Tfreezing 0.0E1				/* freezing temperature [Celsius]*/
#define ka 0.41						/* Von Karman constant*/
#define mu_l 0.001787				/* Dynamic viscosity of water at 0 degrees Celsius*/
#define Asurr 0.0		/* Albedo of sorrounding terrain */
#define wsn_vis 0.8					//snow on canopy: scattering parameters
#define wsn_nir 0.4
#define Bsnd_vis 0.5
#define Bsnd_nir 0.5
#define Bsnb_vis 0.5
#define Bsnb_nir 0.5
#define Tol_h_mount 2.0
#define Tol_h_flat 20.0
#define veg_jumping_exp 3.0


//Meteo data
#define iPt 0		/*Precipitation*/
#define iWs iPt+1	/*Wind speed*/
#define iWd iWs+1	/*Wind direction*/
#define iRh iWd+1	/*Relative humidity*/
#define iT iRh+1	/*Air temperature*/
#define iTlr iT+1	/*Lapse rate*/
#define iPs iTlr+1	/*Air Pressure*/
#define iSW iPs+1	/*global shortwave radiation*/
#define iSWb iSW+1	/*direct SW*/
#define iSWd iSWb+1	/*diffuse SW*/
#define itauC iSWd+1/*sky trasmissivity*/
#define iC	 itauC+1/*Cloudiness*/
#define iLWi iC+1	/*incoming longwave*/
#define iSWn iLWi+1	/*net shortwave*/
#define iLWn iSWn+1	/*net longwave*/
#define iH iLWn+1	/*Sensible heat flux*/
#define iLE iH+1	/*Latent heat flux*/
#define nmet iLE+1	/* number of meteo variables */

//Files
/*#define fpar 1						//parameter file
#define fopt fpar+1					//options
#define fspar fopt+1				//soil parameters
#define fmet fspar+1				//meteo
#define fhor fmet+1					//horizon
#define fdem fhor+1					//digital elevation model
#define fdd fdem+1					//drainage directions
#define flu fdd+1					//land use
#define fsoil flu+1					//soil type map
#define fcurv fsoil+1				//curvature (binary)
#define fsky fcurv+1				//sky view factor
#define fslp fsky+1					//slope
#define fnet fslp+1					//channel network
#define fgrad fnet+1				//gradient along drainage directions
#define fasp fgrad+1				//aspect (0 north, then clockwise)
#define farea fasp+1				//area accounting for slope correction
#define fdist farea+1				//distance from outlet
#define ftca fdist+1				//total contributing areas (in pixels)
#define fsn0 ftca+1					//initial snow depth (mm)
#define fsnag0 fsn0+1				//initial snow age (days)
#define fgl0 fsnag0+1				//initial glacier depth (mm)
#define ferr fgl0+1					//error file
#define fQ ferr+1					//(o.) output discharge file
#define fbas fQ+1					//o. basin variables
#define farank fbas+1				//o. altimetric ranks variables
#define fpoint farank+1				//o. point variables
#define fTz fpoint+1				//o. temperature profiles
#define fpsiz fTz+1					//o. psi profiles
#define fliqz fpsiz+1				//o. water content profiles
#define ficez fliqz+1				//o. ice content profiles
#define fsnz ficez+1				//o. snow data
#define fglz fsnz+1					//o. glacier data
#define fT fglz+1					//o. temperature maps
#define fliq fT+1					//o. water content maps
#define fice fliq+1					//o. ice content maps
#define fhsup fice+1				//o. water over the surface (mm) maps
#define falb fhsup+1				//o. albedo maps
#define fRn falb+1					//o. radiation maps
#define fG fRn+1					//o. surface heat flux maps
#define fH fG+1						//o. sensible heat flux maps
#define fLE fH+1					//o. latent heat flux maps
#define fTs fLE+1					//o. surface temperature maps
#define fprec fTs+1					//o. precipitation maps
#define fcint fprec+1				//o. precipitation intercepted by canopy maps
#define fpsi fcint+1				//o. psi maps
#define fsn fpsi+1					//o. snow maps
#define fgl fsn+1					//o. glacier maps
#define fmsn fgl+1					//o. snow melted maps
#define fssn fmsn+1					//o. snow sublimated maps
#define fmgl fssn+1					//o. glacier ice melted maps
#define fsgl fmgl+1					//o. glacier ice sublimated maps
#define fSW fsgl+1					//o. shortwave radiation maps
#define fTa fSW+1					//o. air temperature maps
#define fwspd fTa+1					//o. wind speed maps
#define fwdir fwspd+1				//o. wind direction maps
#define frh fwdir+1					//o. relative humidity maps
#define fsnd frh+1					//o. snow density maps
#define fgld fsnd+1					//o. glacier ice density maps
#define fsndur fgld+1				//o. snow duration maps (days)
#define fsnav fsndur+1				//o. averaged snow depth
#define fmeltlu fsnav+1				//o. precipitation/snow melted/glacier melted for each land use type
#define rpsi fmeltlu+1					//recover file (f.) psi
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
#define rns rTi+1					//r. number of snow layers
#define rni rns+1					//r. number of glacier layers
#define rsnag rni+1					//r. snow age
#define rhsup rsnag+1				//r. water over the surface
#define rQch rhsup+1					//r. water stored in channels
#define nfiles rQch					//number of files
*/
//soil data
#define jdz 1						//layer thickness [mm]
#define jpsi jdz+1					//initial psi [mm]
#define jT jpsi+1					//initial temperature [C]
#define jKh jT+1					//lateral hydr. conductivity [mm/s]
#define jKv jKh+1					//vertical hydr. conductivity [mm/s]
#define jres jKv+1					//residual wat.cont.
#define jsat jres+1					//porosity
#define ja jsat+1					//alpha[mm^-1]
#define jns ja+1					//n
#define jv jns+1					//v
#define jkt jv+1					//thermal conductivity
#define jct jkt+1					//thermal capacity
#define jpsimin jct+1				//minimum psi (for Richards' eq. integration purposes)
#define jlatfl jpsimin+1			//flag to allow lateral flow
#define jsf jlatfl+1				//flag to allow soil freezing
#define jKav jsf+1					//flag for which type of average hydr.conductivity (harmonic or arithmetic)
#define nsoilprop jKav				//number of soil properties considered

//land use data
#define jz0soil 1					//roughness length for soil
#define jz0thressoil jz0soil+1		//threshold on snow depth to change roughness length to snow covered values in soil area
#define jHveg jz0thressoil+1		//vegetation height
#define jz0veg jHveg+1				//roughness length for vegetation
#define jd0 jz0veg+1				//displacement height for vegetation
#define jz0thresveg jd0+1			//threshold on snow depth to change roughness length to snow covered values in vegetated area
#define jLAIw jz0thresveg+1			//LAI winter
#define jLAIs jLAIw+1				//LAI summer
#define jcf jLAIs+1					//Canopy gap fraction (k)
#define jroot jcf+1					//root depth [mm]
#define jrs jroot+1					//canopy transpiration coefficient
#define jtwp jrs+1					//wilting point water cont.
#define jtfc jtwp+1					//field capacity water cont.
#define jvholdsn jtfc+1				//vegetation holding snow capacity (used for blowing snow)
#define jvR_vis jvholdsn+1			//vegetation in the visible spectrum
#define jvR_nir jvR_vis+1			//vegetation in the near infrared spectrum
#define jvT_vis jvR_nir+1			//vegetation in the visible spectrum
#define jvT_nir jvT_vis+1			//vegetation in the near infrared spectrum
#define jvCh jvT_nir+1				//departure of leaf angles from a random distribution (1 horizontal, 0 random, -1 vertical)
#define jcd jvCh+1					//surface density of canopy [kg/(m2*LAI)]
#define ja_vis jcd+1				//ground albedo in the visible spectrum
#define ja_nir ja_vis+1				//ground albedo in the near infrared spectrum
#define jemg ja_nir+1				//soil emissivity
#define jcm jemg+1					//gauckler strickler (1/manning) coefficient
#define jzb jcm+1					//Depth 0 Temperature amplitude
#define jtb jzb+1					//Temperature at the Depth above
#define nlandprop jtb				//number of land use properties
