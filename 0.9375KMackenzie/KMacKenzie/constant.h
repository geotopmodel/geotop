
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion MacLavagna

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 MacLavagna.
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
#define c_liq 4188.0				/* heat capacity of water		[J/(kg/K)]*/
#define c_ice 2117.0				/* heat capacity of ice		[J/(kg/K)]*/
#define CA 0.34						/* tunable parameter in order to have surface temperature instead layer center temperature*/
#define KNe 0.50					/* Krank-Nicholson parameter for egy balance*/
#define KNw 1.00					/* Krank-Nicholson parameter for water balance*/
#define Tfreezing 0.0E1				/* freezing temperature [Celsius]*/
#define ka 0.41						/* Von Karman constant*/
#define mu_l 0.001787				/* Dynamic viscosity of water at 0 degrees Celsius*/
#define Asurr 0.8

//Meteo data
#define iPt 0						/*Precipitation*/
#define iWs iPt+1					/*Wind speed*/
#define iWd iWs+1					/*Wind direction*/
#define iRh iWd+1					/*Relative humidity*/
#define iT iRh+1					/*Air temperature*/
#define iTlr iT+1					/*Lapse rate*/
#define iPs iTlr+1					/*Air Pressure*/
#define iSW iPs+1					/*global shortwave radiation*/
#define iSWb iSW+1					/*direct SW*/
#define iSWd iSWb+1					/*diffuse SW*/
#define iC	 iSWd+1					/*Cloudiness*/
#define iSWi iC+1					/*incoming shortwave*/
#define iLWi iSWi+1					/*incoming longwave*/
#define iSWo iLWi+1					/*outgoing shortwave*/
#define iLWo iSWo+1					/*outgoing longwave*/
#define iH iLWo+1					/*Sensible heat flux*/
#define iLE iH+1					/*Latent heat flux*/
#define nmet iLE+1

/* to removed
 ///Files
#define fpar 1
#define fopt fpar+1
#define fspar fopt+1
#define fmet fspar+1
#define fhor fmet+1
#define fdem fhor+1
#define fdd fdem+1
#define flu fdd+1
#define fsoil flu+1
#define fcurv fsoil+1
#define fsky fcurv+1
#define fslp fsky+1
#define fnet fslp+1
#define fgrad fnet+1
#define fasp fgrad+1
#define farea fasp+1
#define fdist farea+1
#define ftca fdist+1
#define fsn0 ftca+1
#define fsnag0 fsn0+1
#define fgl0 fsnag0+1
#define ferr fgl0+1
#define fQ ferr+1
#define fbas fQ+1
#define farank fbas+1
#define fpoint farank+1
#define fTz fpoint+1
#define fpsiz fTz+1
#define fliqz fpsiz+1
#define ficez fliqz+1
#define fsnz ficez+1
#define fglz fsnz+1
#define fT fglz+1
#define fliq fT+1
#define fice fliq+1
#define fhsup fice+1
#define falb fhsup+1
#define fRn falb+1
#define fG fRn+1
#define fH fG+1
#define fLE fH+1
#define fTs fLE+1
#define fprec fTs+1
#define fcint fprec+1
#define fpsi fcint+1
#define fsn fpsi+1
#define fgl fsn+1
#define fmsn fgl+1
#define fssn fmsn+1
#define fmgl fssn+1
#define fsgl fmgl+1
#define fSW fsgl+1
#define fTa fSW+1
#define fwspd fTa+1
#define fwdir fwspd+1
#define frh fwdir+1
#define fsnd frh+1
#define fgld fsnd+1
#define fsndur fgld+1
#define fsnav fsndur+1
#define fmeltlu fsnav+1
#define pH fmeltlu+1
#define pLE pH+1
#define pSWin pLE+1
#define pSWout pSWin+1
#define pLWin pSWout+1
#define pLWout pLWin+1
#define pTs pLWout+1
#define pTa pTs+1
#define pVspd pTa+1
#define pVdir pVspd+1
#define pRH pVdir+1
#define pD pRH+1
#define pth pD+1
#define rpsi pth+1
#define riceg rpsi+1
#define rTg riceg+1
#define rDzs rTg+1
#define rwls rDzs+1
#define rwis rwls+1
#define rTs rwis+1
#define rDzi rTs+1
#define rwli rDzi+1
#define rwii rwli+1
#define rTi rwii+1
#define rns rTi+1
#define rni rns+1
#define rsnag rni+1
#define rhsup rsnag+1
#define rwt rhsup+1
#define rQch rwt+1
#define rSFA rQch+1
#define fHpatch rSFA+1
#define nfiles fHpatch

 */
//soil data
#define jdz 1
#define jpsi jdz+1
#define jT jpsi+1
#define jKh jT+1
#define jKv jKh+1
#define jres jKv+1
#define jsat jres+1
#define ja jsat+1
#define jns ja+1
#define jv jns+1
#define jkt jv+1
#define jct jkt+1
#define jpsimin jct+1
#define jlatfl jpsimin+1
#define jsf jlatfl+1
#define jKav jsf+1
#define nsoilprop jKav

//land use data
#define jz0 1
#define jz0zt jz0+1
#define jd0 jz0zt+1
#define jz0thres jd0+1
#define jhc jz0thres+1
#define jfc jhc+1
#define jLAIw jfc+1
#define jLAIs jLAIw+1
#define jroot jLAIs+1
#define jrs jroot+1
#define jtwp jrs+1
#define jtfc jtwp+1
#define jvholdsn jtfc+1
#define jalbedo jvholdsn+1
#define jem jalbedo+1
#define jcm jem+1
#define jGbottom jcm+1
#define nlandprop jGbottom
