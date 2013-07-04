
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


//Missing file string
#define LOCAL_MISSING_FILE "MISSING_FILE"

//Files
#define fcontrol 1					//control parameters
#define fpar fcontrol+1				//parameter file
#define fopt fpar+1					//options
#define fspar fopt+1				//soil parameters
#define fspar2 fspar+1				//soil bedrock parameter
#define fmet fspar2+1				//meteo
#define fLRs fmet+1					//lapse rates
#define fhor fLRs+1					//horizon
#define fvegpar fhor+1				//vegetation parameter
#define fdem fvegpar+1				//digital elevation model (m)
#define fdd fdem+1					//drainage directions
#define flu fdd+1					//land use
#define fsoil flu+1					//soil type map
#define fsky fsoil+1				//sky view factor
#define fslp fsky+1					//slope
#define fnet fslp+1					//channel network
#define fasp fnet+1					//aspect (0 north, then clockwise)
#define fbed fasp+1					//bedrock topography (m)
#define fsn0 fbed+1					//initial snow depth (mm)
#define fsnag0 fsn0+1				//initial snow age (days)
#define fgl0 fsnag0+1				//initial glacier depth (mm)
#define ferr fgl0+1					//error file
#define fQ ferr+1					//(o.) output discharge file
#define fbas fQ+1					//o. basin variables
#define fpoint fbas+1				//o. point variables
#define fTz fpoint+1				//o. temperature profiles
#define fpsiz fTz+1					//o. psi profiles
#define fpsiztot fpsiz+1			//o. psi profiles
#define fliqz fpsiztot+1			//o. water content profiles
#define ficez fliqz+1				//o. ice content profiles
#define fsnz ficez+1				//o. snow data
#define fglz fsnz+1					//o. glacier data
#define fSCA fglz+1					//file giving the fraction of snow free areas and corresponding properties 
#define fT fSCA+1					//o. temperature maps
#define fliq fT+1					//o. water content maps
#define fice fliq+1					//o. ice content maps
#define fliqsup fice+1				//o. water content in the soil at the surface
#define fhsup fliqsup+1				//o. water over the surface (mm) maps
#define fRn fhsup+1					//o. radiation maps
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
#define fswe frh+1					//o. snow density maps
#define fgwe fswe+1					//o. glacier ice density maps
#define fsndur fgwe+1				//o. snow duration maps (hrs)
#define fsnav fsndur+1				//o. averaged snow depth 
#define fthawed fsnav+1				//thawed soil depth
#define fwtable fthawed+1			//water table depth
#define fftable fwtable+1			//frost table depth
#define pH fftable+1
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
#define rns rTi+1					//r. number of snow layers
#define rni rns+1					//r. number of glacier layers
#define rsnag_adim rni+1			//r. snow age
#define rsnag_dim rsnag_adim+1
#define rhsup rsnag_dim+1			//r. water over the surface 
#define rwcrn rhsup+1				//r. water stored on canopy
#define rwcsn rwcrn+1
#define rTv rwcsn+1
#define rQch rTv+1					//r. water stored in channels
#define nfiles rQch					//number of files
