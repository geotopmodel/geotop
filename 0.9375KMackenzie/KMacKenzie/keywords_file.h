
/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion KMackenzie

Copyright, 2008 Stefano Endrizzi, Emanuele Cordano, Riccardo Rigon, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 KMackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOtop is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/




#define I_CONTROL_PARAMETERS   1        /* text files cotaiig parameters releted to the simulation (timestep, format output, etc previously contained in the geotop.inpts file of the old versions) (without extension) */
#define	I_PARAMETERS	     I_CONTROL_PARAMETERS+1	     /*	  text file containing parameters necessary for simulations (without extension)	*/
#define	I_OPTIONS	         I_PARAMETERS+1	/*	  text file containing run options (without extension)	*/
#define	I_SOIL_INFO	         I_OPTIONS+1	/*	 text file containing soil properties for each soil type (without extension)	*/
#define	I_METEO_PRIMARY	     I_SOIL_INFO+1	/*	 text file containing time series of meteo data for each meteo station (without extension and meteo station number)	*/
#define	I_METEO_HORIZON	     I_METEO_PRIMARY+1	/*	 text file containg the horizons for each stations (without extension and meteo station number)	*/
#define I_METEO_CLOUDCOVER   I_METEO_HORIZON+1    /* text file containing cloud cover informatio (old "ii_cloud.txt")(without extension) */
#define	I_MORPHO_ELEVATION	 I_METEO_CLOUDCOVER+1 	/*	" 		Digital Terren Model of the basin  (without extension)"	*/
#define	I_MORPHO_DD	         I_MORPHO_ELEVATION+1	/*	" 				Map of  Drainage Diractions (without extension) "	*/
#define	I_LANDUSE_MAP	     I_MORPHO_DD+1	/*	 map of land use (without extension) 	*/
#define	I_SOILTYPE_MAP	     I_LANDUSE_MAP+1	/*	 map of soil type (without extension) 	*/
#define	I_MORPHO_NABLA	           I_SOILTYPE_MAP+1	/*	"  			Map of Curvature (Nabla) (0 CONCAVE ZONES - 1 CONVEX ZONES))  (without extension)  "	*/
#define	I_MORPHO_SKYVIEWFACTOR	   I_MORPHO_NABLA+1	/*	     Map of Sky View Factor (without extension)  	*/
#define	I_MORPHO_SLOPE	           I_MORPHO_SKYVIEWFACTOR+1	/*	"  			Map of Slope (gradient) along Maximum Slope Direction (without extension) "	*/
#define	I_MORPHO_CHANNELNETWORK	   I_MORPHO_SLOPE+1	/*	" 	Map of Channel Network  (without extension) "	*/
#define	I_MORPHO_DD_SLOPE	       I_MORPHO_CHANNELNETWORK+1	/*	" 			Map of Slope along Drainage Directions   (without extension)"	*/
#define	I_MORPHO_ASPECT	           I_MORPHO_DD_SLOPE+1	/*	"  			Map of Aspect  (without extension)"	*/
#define	I_MORPHO_AREA	           I_MORPHO_ASPECT+1	/*	"   			Map of Total Contributing Area taking into account the slope [square meters] (without extension) "	*/
#define	I_MORPHO_PIXELDISTANCE	   I_MORPHO_AREA+1	/*	     Map of Distance from pixels to outlet (without extension)	*/
#define	I_MORPHO_TCA	           I_MORPHO_PIXELDISTANCE+1	/*	"   			Map of Total Contributing Area expressed in number of pixels (without extension) "	*/
#define	I_CONDINIT_SNOWDEPTH	   I_MORPHO_TCA+1	/*	       map of snow depth (without extension)	*/
#define	I_SNOWPROP_AGE	           I_CONDINIT_SNOWDEPTH+1	/*	             map of snow age (jn days)  (without extension)	*/
#define	I_CONDINIT_GLACIERDEPTH	   I_SNOWPROP_AGE+1	/*	    map of glacier depth (without extension)	*/
#define	O_ERRORS	               I_CONDINIT_GLACIERDEPTH+1	/*	                   output file of errors (without extension)	*/
#define	O_WATERDISCHARGE_AT_OUTLET O_ERRORS+1	/*	 output text file with water discharge at the outlet (without extension)	*/
#define	O_TIMESERIES_BASIN	              O_WATERDISCHARGE_AT_OUTLET+1	/*	         output text file with time series of the mean variables of the whole basin (whitout extention) 	*/
#define	O_TIMESERIES_ES	                  O_TIMESERIES_BASIN+1	/*	            output text files with time series of the mean variables for each altimetric stripes (whitout extention) 	*/
#define	O_TIMESERIES_CONTROLPIXELS	      O_TIMESERIES_ES+1	/*	 output text files with time series of the mean variables at the control pixels (whitout extention) 	*/
#define	O_SOILTEMPERTURE_CONTROLPIXELS	  O_TIMESERIES_CONTROLPIXELS+1	/*	 output text files with time series of soil temperature at the control pixels (whitout extantion) 	*/
#define	O_PRESSUREHEAD_CONTROLPIXELS	  O_SOILTEMPERTURE_CONTROLPIXELS+1	/*	 output text files with time series of soil water pressure head at the control pixels (whitout extention) 	*/
#define	O_WATERCONTENT_CONTROLPIXELS	  O_PRESSUREHEAD_CONTROLPIXELS+1	/*	 output text files with time series of soil water content at the control pixels (whitout extention) 	*/
#define	O_ICECONTENT_CONTROLPIXELS	      O_WATERCONTENT_CONTROLPIXELS+1	/*	 output text files with time series of volumetric ice content at the control pixels (whitout extention) 	*/
#define	O_SNOWPROFILE_CONTROLPIXELS	      O_ICECONTENT_CONTROLPIXELS+1	/*	  output text files with time series of snow profiles at the control pixels (whitout extention) 	*/
#define	O_GLACIERPROFILE_CONTROLPIXELS	  O_SNOWPROFILE_CONTROLPIXELS+1	/*	  output text files with time series of glacier profiles at the control pixels (whitout extention) 	*/
#define	O_SOILTEMPERATURE_DISTRIBUTED_TENSOR	 O_GLACIERPROFILE_CONTROLPIXELS+1	/*	  output 3D plus time tensor for soil tempeture distributed  in the whole basin (whitout extention) 	*/
#define	O_WATERCONTENT_DISTRIBUTED_TENSOR	     O_SOILTEMPERATURE_DISTRIBUTED_TENSOR+1	/*	     output 3D plus time tensor for soil water content  distributed  in the whole basin (whitout extention) 	*/
#define	O_ICECONTENT_DISTRIBUTED_TENSOR	         O_WATERCONTENT_DISTRIBUTED_TENSOR+1	/*	       output 3D plus time tensor for soil ice content  distributed  in the whole basin (whitout extention) 	*/
#define	O_SURFACE_WATER_DEPTH_MAP	             O_ICECONTENT_DISTRIBUTED_TENSOR+1	/*	 output 2D plus time tensor for water surface depth (h_sup) [mm] distributed  in the whole basin (whitout extention) 	*/
#define	O_ALBEDO_WITH_SNOW_DISTRIBUTED_MAP	     O_SURFACE_WATER_DEPTH_MAP+1	/*	 output 2D plus time tensor for albedo (taking into account snow distribution)  distributed  in the whole basin (whitout extention) 	*/
#define	O_NET_RADIATION_DISTRIBUTED_MAP	         O_ALBEDO_WITH_SNOW_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for net radiation  distributed  in the whole basin (whitout extention) 	*/
#define	O_GROUND_HEAT_FLUX_DISTRIBUTED_MAP	     O_NET_RADIATION_DISTRIBUTED_MAP+1	/*	  output 2D plus time tensor for ground heat flux distributed  in the whole basin (whitout extention) 	*/
#define	O_SENSIBLE_HEAT_FLUX_DISTRIBUTED_MAP	 O_GROUND_HEAT_FLUX_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for sensible heat flux distributed  in the whole basin (whitout extention) 	*/
#define	O_LATENT_HEAT_FLUX_DISTRIBUTED_MAP	     O_SENSIBLE_HEAT_FLUX_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for latent heat flux distributed  in the whole basin (whitout extention) 	*/
#define	O_SURFACE_TEMEPERATURE_DISTRIBUTED_MAP	 O_LATENT_HEAT_FLUX_DISTRIBUTED_MAP+1	/*	  output 2D plus time tensor for surface temperature distributed  in the whole basin (whitout extention) 	*/
#define	O_PRECIPITATION_DISTRIBUTED_MAP	         O_SURFACE_TEMEPERATURE_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for precipitation distributed  in the whole basin (whitout extention) 	*/
#define	O_INTERCEPTATION_DISTRIBUTED_MAP	O_PRECIPITATION_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for interceptation distributed  in the whole basin (whitout extention) 	*/
#define	O_PRESSUREHEAD_DISTRIBUTED_TENSOR	O_INTERCEPTATION_DISTRIBUTED_MAP+1	/*	     output 3D plus time tensor for soil water pressure head distributed  in the whole basin (whitout extention) 	*/
#define	O_SNOW_WATER_EQUIVALENT_DISTRIBUTED_MAP	O_PRESSUREHEAD_DISTRIBUTED_TENSOR+1	/*	  output 2D plus time tensor for snow water equivalent distributed  in the whole basin (whitout extention) 	*/
#define	O_GLACIER_WATER_EQUIVALENT_DISTRIBUTED_MAP	O_SNOW_WATER_EQUIVALENT_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for glacier water equivalent distributed  in the whole basin (whitout extention) 	*/
#define	O_MELTED_SNOW_DISTRIBUTED_MAP	            O_GLACIER_WATER_EQUIVALENT_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for melted snow water equivalent distributed in the whole basin (whitout extention) 	*/
#define	O_SUBLIMATED_SNOW_DISTRIBUTED_MAP	        O_MELTED_SNOW_DISTRIBUTED_MAP+1	/*	  output 2D plus time tensor for sublimated snow distributed in the whole basin (whitout extention) 	*/
#define	O_MELTED_GLACIER_DISTRIBUTED_MAP	        O_SUBLIMATED_SNOW_DISTRIBUTED_MAP+1	/*	   output 2D plus time tensor for melted ice distributed  in the whole basin (whitout extention) 	*/
#define	O_SUBLIMATED_GLACIER_DISTRIBUTED_MAP	    O_MELTED_GLACIER_DISTRIBUTED_MAP+1	/*	   output 2D plus time tensor for sublimated ice distributed  in the whole basin (whitout extention) 	*/
#define	O_INCOMING_SW_RADIATION_DISTRIBUTED_MAP	    O_SUBLIMATED_GLACIER_DISTRIBUTED_MAP+1	/*	   output 2D plus time tensor for incoming SW radiation distributed  in the whole basin (whitout extention) 	*/
#define	O_AIR_TEMPERATURE_DISTRIBUTED_MAP	        O_INCOMING_SW_RADIATION_DISTRIBUTED_MAP+1	/*	  output 2D plus time tensor for air temperature distributed  in the whole basin (whitout extention) 	*/
#define	O_WINDVELOCTY_DISTRIBUTED_MAP	            O_AIR_TEMPERATURE_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for for Wind Velocity distributed  in the whole basin (whitout extension) 	*/
#define	O_WINDDIRECTION_DISTRIBUTED_MAP	            O_WINDVELOCTY_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for for Wind Directions distributed  in the whole basin (whitout extension) 	*/
#define	O_RELHUMIDITY_DISTRIBUTED_MAP	            O_WINDDIRECTION_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for for Relaive Humidity distributed  in the whole basin (whitout extension) 	*/
#define	O_SNOWDENSITY_DISTRIBUTED_MAP	            O_RELHUMIDITY_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for for Snow Density distributed  in the whole basin (whitout extension) 	*/
#define	O_GLACDENSITY_DISTRIBUTED_MAP	            O_SNOWDENSITY_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for for Glacier Density distributed  in the whole basin (whitout extension) 	*/
#define	O_SNOW_DURATION_DISTRIBUTED_MAP	            O_GLACDENSITY_DISTRIBUTED_MAP+1	/*	  output 2D plus time tensor for snow duration distributed  in the whole basin (whitout extention) 	*/
#define	O_AVERAGE_SNOW_FROM_BEGGIN_DISTRIBUTED_MAP	O_SNOW_DURATION_DISTRIBUTED_MAP+1	/*	 output 2D plus time tensor for averaged snow since begin of simulation distributed  in the whole basin (whitout extention) 	*/
#define	O_MELT_FLUXES_vs_LANDUSE	                O_AVERAGE_SNOW_FROM_BEGGIN_DISTRIBUTED_MAP+1	/*	  output textfile for melt fluxes for each lnd use class the whole basin (whitout extention) 	*/
#define	O_K_pH	                                    O_MELT_FLUXES_vs_LANDUSE+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pLE	                                    O_K_pH+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pSWin	                                O_K_pLE+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pSWout	                                O_K_pSWin+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pLWin	                                O_K_pSWout+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pLWout	                                O_K_pLWin+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pTs	                                    O_K_pLWout+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pTa	                                    O_K_pTs+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pVspd	                                O_K_pTa+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pVdir	                                O_K_pVspd+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pRH	                                    O_K_pVdir+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_pD	                                    O_K_pRH+1	/*	 to be commented by Stefano Endrizzi	*/
#define	O_K_ptheta	                                O_K_pD+1	/*	 to be commented by Stefano Endrizzi	*/
#define	REC_PRESSUREHEAD	                        O_K_ptheta+1	/*	   recover file for pressure head (without extention)	*/
#define	REC_ICECONTENT	                            REC_PRESSUREHEAD+1	/*	      recover file for ice content (without extention)	*/
#define	REC_TEMPERATURE	                            REC_ICECONTENT+1	/*	     recover file for temparature (without extention)	*/
#define	REC_SNOWLAYER	                            REC_TEMPERATURE+1	/*	       recover file for snow layer (without extention)	*/
#define	REC_WLIQ_SNOW	                            REC_SNOWLAYER+1	/*	       recover file for melted liquid water content WLIQ (without extention)	*/
#define	REC_WICE_SNOW	                            REC_WLIQ_SNOW+1	/*	       recover file for melted solid water content WICE (without extention)	*/
#define	REC_SNOW_TEMPERATURE	                    REC_WICE_SNOW+1	/*	 recover file for snow temperature (without extention)	*/
#define	REC_GLACIERLAYER	                        REC_SNOW_TEMPERATURE+1	/*	    recover file for glacier layer (without extention)	*/
#define	REC_WLIQ_GLACIER	                        REC_GLACIERLAYER+1	/*	   recover file for WLIQ in glacier (without extention)	*/
#define	REC_WICE_GLACIER	                        REC_WLIQ_GLACIER+1	/*	   recover file for WICE in glacier (without extention)	*/
#define	REC_GLACIER_TEMPERATURE	                    REC_WICE_GLACIER+1	/*	   recover file for glacier temperature (without extention)	*/
#define	REC_NUMBER_SNOWLAYER	                    REC_GLACIER_TEMPERATURE+1	/*	      recover file for number of snow layers(without extention)	*/
#define	REC_NUMBER_GLACIERLAYER	                    REC_NUMBER_SNOWLAYER+1	/*	     recover file for number of glacier layers(without extention)	*/
#define	REC_DIMENSIONLESS_SNOW_AGE	                REC_NUMBER_GLACIERLAYER+1	/*	  recover file for dimensionless snow age (without extention)	*/
#define	REC__h_sup	                                REC_DIMENSIONLESS_SNOW_AGE+1	/*	 recover file for h_sup (without extention)	*/
#define	REC__wt	                                    REC__h_sup+1	/*	 recover file for wt (without extention)	*/
#define	REC__Qchannel	                            REC__wt+1	/*	 recover file for Qchannel (without extention)	*/
#define	REC_SFA	                                    REC__Qchannel+1	/*	 recover file for SFR (without extention)	*/
#define	zSFA	                                    REC_SFA+1	/*	 recover file for zSFR (without extention)	*/
#define	nfiles  									zSFA /*number of files */


/* alias Stefano */
#define	fpar	I_PARAMETERS	/*	 text file containing parameters necessary for simulations (without extension)	*/
#define	fopt	I_OPTIONS	/*	  text file containing run options (without extension)	*/
#define	fspar	I_SOIL_INFO	/*	 text file containing soil properties for each soil type (without extension)	*/
#define	fmet	I_METEO_PRIMARY	/*	 text file containing time series of meteo data for each meteo station (without extension and meteo station number)	*/
#define	fhor	I_METEO_HORIZON	/*	 text file containg the horizons for each stations (without extension and meteo station number)	*/
#define	fdem	I_MORPHO_ELEVATION	/*	" 		Digital Terren Model of the basin  (without extension)"	*/
#define	fdd	I_MORPHO_DD	/*	" 				Map of  Drainage Diractions (without extension) "	*/
#define	flu	I_LANDUSE_MAP	/*	 map of land use (without extension) 	*/
#define	fsoil	I_SOILTYPE_MAP	/*	 map of soil type (without extension) 	*/
#define	fcurv	I_MORPHO_NABLA	/*	"  			Map of Curvature (Nabla) (0 CONCAVE ZONES - 1 CONVEX ZONES))  (without extension)  "	*/
#define	fsky	I_MORPHO_SKYVIEWFACTOR	/*	     Map of Sky View Factor (without extension)  	*/
#define	fslp	I_MORPHO_SLOPE	/*	"  			Map of Slope (gradient) along Maximum Slope Direction (without extension) "	*/
#define	fnet	I_MORPHO_CHANNELNETWORK	/*	" 	Map of Channel Network  (without extension) "	*/
#define	fgrad	I_MORPHO_DD_SLOPE	/*	" 			Map of Slope along Drainage Directions   (without extension)"	*/
#define	fasp	I_MORPHO_ASPECT	/*	"  			Map of Aspect  (without extension)"	*/
#define	farea	I_MORPHO_AREA	/*	"   			Map of Total Contributing Area taking into account the slope [square meters] (without extension) "	*/
#define	fdist	I_MORPHO_PIXELDISTANCE	/*	     Map of Distance from pixels to outlet (without extension)	*/
#define	ftca	I_MORPHO_TCA	/*	"   			Map of Total Contributing Area expressed in number of pixels (without extension) "	*/
#define	fsn0	I_CONDINIT_SNOWDEPTH	/*	       map of snow depth (without extension)	*/
#define	fsnag0	I_SNOWPROP_AGE	/*	             map of snow age (jn days)  (without extension)	*/
#define	fgl0	I_CONDINIT_GLACIERDEPTH	/*	    map of glacier depth (without extension)	*/
#define	ferr	O_ERRORS	/*	                   output file of errors (without extension)	*/
#define	fQ	O_WATERDISCHARGE_AT_OUTLET	/*	 output text file with water discharge at the outlet (without extension)	*/
#define	fbas	O_TIMESERIES_BASIN	/*	         output text file with time series of the mean variables of the whole basin (whitout extention) 	*/
#define	farank	O_TIMESERIES_ES	/*	            output text files with time series of the mean variables for each altimetric stripes (whitout extention) 	*/
#define	fpoint	O_TIMESERIES_CONTROLPIXELS	/*	 output text files with time series of the mean variables at the control pixels (whitout extention) 	*/
#define	fTz	O_SOILTEMPERTURE_CONTROLPIXELS	/*	 output text files with time series of soil temperature at the control pixels (whitout extantion) 	*/
#define	fpsiz	O_PRESSUREHEAD_CONTROLPIXELS	/*	 output text files with time series of soil water pressure head at the control pixels (whitout extention) 	*/
#define	fliqz	O_WATERCONTENT_CONTROLPIXELS	/*	 output text files with time series of soil water content at the control pixels (whitout extention) 	*/
#define	ficez	O_ICECONTENT_CONTROLPIXELS	/*	 output text files with time series of volumetric ice content at the control pixels (whitout extention) 	*/
#define	fsnz	O_SNOWPROFILE_CONTROLPIXELS	/*	  output text files with time series of snow profiles at the control pixels (whitout extention) 	*/
#define	fglz	O_GLACIERPROFILE_CONTROLPIXELS	/*	  output text files with time series of glacier profiles at the control pixels (whitout extention) 	*/
#define	fT	O_SOILTEMPERATURE_DISTRIBUTED_TENSOR	/*	  output 3D plus time tensor for soil tempeture distributed  in the whole basin (whitout extention) 	*/
#define	fliq	O_WATERCONTENT_DISTRIBUTED_TENSOR	/*	     output 3D plus time tensor for soil water content  distributed  in the whole basin (whitout extention) 	*/
#define	fice	O_ICECONTENT_DISTRIBUTED_TENSOR	/*	       output 3D plus time tensor for soil ice content  distributed  in the whole basin (whitout extention) 	*/
#define	fhsup	O_SURFACE_WATER_DEPTH_MAP	/*	 output 2D plus time tensor for water surface depth (h_sup) [mm] distributed  in the whole basin (whitout extention) 	*/
#define	falb	O_ALBEDO_WITH_SNOW_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for albedo (taking into account snow distribution)  distributed  in the whole basin (whitout extention) 	*/
#define	fRn	O_NET_RADIATION_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for net radiation  distributed  in the whole basin (whitout extention) 	*/
#define	fG	O_GROUND_HEAT_FLUX_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for ground heat flux distributed  in the whole basin (whitout extention) 	*/
#define	fH	O_SENSIBLE_HEAT_FLUX_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for sensible heat flux distributed  in the whole basin (whitout extention) 	*/
#define	fLE	O_LATENT_HEAT_FLUX_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for latent heat flux distributed  in the whole basin (whitout extention) 	*/
#define	fTs	O_SURFACE_TEMEPERATURE_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for surface temperature distributed  in the whole basin (whitout extention) 	*/
#define	fprec	O_PRECIPITATION_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for precipitation distributed  in the whole basin (whitout extention) 	*/
#define	fcint	O_INTERCEPTATION_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for interceptation distributed  in the whole basin (whitout extention) 	*/
#define	fpsi	O_PRESSUREHEAD_DISTRIBUTED_TENSOR	/*	     output 3D plus time tensor for soil water pressure head distributed  in the whole basin (whitout extention) 	*/
#define	fsn	O_SNOW_WATER_EQUIVALENT_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for snow water equivalent distributed  in the whole basin (whitout extention) 	*/
#define	fgl	O_GLACIER_WATER_EQUIVALENT_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for glacier water equivalent distributed  in the whole basin (whitout extention) 	*/
#define	fmsn	O_MELTED_SNOW_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for melted snow water equivalent distributed in the whole basin (whitout extention) 	*/
#define	fssn	O_SUBLIMATED_SNOW_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for sublimated snow distributed in the whole basin (whitout extention) 	*/
#define	fmgl	O_MELTED_GLACIER_DISTRIBUTED_MAP	/*	   output 2D plus time tensor for melted ice distributed  in the whole basin (whitout extention) 	*/
#define	fsgl	O_SUBLIMATED_GLACIER_DISTRIBUTED_MAP	/*	   output 2D plus time tensor for sublimated ice distributed  in the whole basin (whitout extention) 	*/
#define	fSW	O_INCOMING_SW_RADIATION_DISTRIBUTED_MAP	/*	   output 2D plus time tensor for incoming SW radiation distributed  in the whole basin (whitout extention) 	*/
#define	fTa	O_AIR_TEMPERATURE_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for air temperature distributed  in the whole basin (whitout extention) 	*/
#define	fwspd	O_WINDVELOCTY_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for for Wind Velocity distributed  in the whole basin (whitout extension) 	*/
#define	fwdir	O_WINDDIRECTION_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for for Wind Directions distributed  in the whole basin (whitout extension) 	*/
#define	frh	O_RELHUMIDITY_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for for Relaive Humidity distributed  in the whole basin (whitout extension) 	*/
#define	fsnd	O_SNOWDENSITY_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for for Snow Density distributed  in the whole basin (whitout extension) 	*/
#define	fgld	O_GLACDENSITY_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for for Glacier Density distributed  in the whole basin (whitout extension) 	*/
#define	fsndur	O_SNOW_DURATION_DISTRIBUTED_MAP	/*	  output 2D plus time tensor for snow duration distributed  in the whole basin (whitout extention) 	*/
#define	fsnav	O_AVERAGE_SNOW_FROM_BEGGIN_DISTRIBUTED_MAP	/*	 output 2D plus time tensor for averaged snow since begin of simulation distributed  in the whole basin (whitout extention) 	*/
#define	fmeltlu	O_MELT_FLUXES_vs_LANDUSE	/*	  output textfile for melt fluxes for each lnd use class the whole basin (whitout extention) 	*/
#define	pH	O_K_pH	/*	 to be commented by Stefano Endrizzi	*/
#define	pLE	O_K_pLE	/*	 to be commented by Stefano Endrizzi	*/
#define	pSWin	O_K_pSWin	/*	 to be commented by Stefano Endrizzi	*/
#define	pSWout	O_K_pSWout	/*	 to be commented by Stefano Endrizzi	*/
#define	pLWin	O_K_pLWin	/*	 to be commented by Stefano Endrizzi	*/
#define	pLWout	O_K_pLWout	/*	 to be commented by Stefano Endrizzi	*/
#define	pTs	O_K_pTs	/*	 to be commented by Stefano Endrizzi	*/
#define	pTa	O_K_pTa	/*	 to be commented by Stefano Endrizzi	*/
#define	pVspd	O_K_pVspd	/*	 to be commented by Stefano Endrizzi	*/
#define	pVdir	O_K_pVdir	/*	 to be commented by Stefano Endrizzi	*/
#define	pRH	O_K_pRH	/*	 to be commented by Stefano Endrizzi	*/
#define	pD	O_K_pD	/*	 to be commented by Stefano Endrizzi	*/
#define	pth	O_K_ptheta	/*	 to be commented by Stefano Endrizzi	*/
#define	rpsi	REC_PRESSUREHEAD	/*	   recover file for pressure head (without extention)	*/
#define	riceg	REC_ICECONTENT	/*	      recover file for ice content (without extention)	*/
#define	rTg	REC_TEMPERATURE	/*	     recover file for temparature (without extention)	*/
#define	rDzs	REC_SNOWLAYER	/*	       recover file for snow layer (without extention)	*/
#define	rwls	REC_WLIQ_SNOW	/*	       recover file for melted liquid water content WLIQ (without extention)	*/
#define	rwis	REC_WICE_SNOW	/*	       recover file for melted solid water content WICE (without extention)	*/
#define	rTs	REC_SNOW_TEMPERATURE	/*	 recover file for snow temperature (without extention)	*/
#define	rDzi	REC_GLACIERLAYER	/*	    recover file for glacier layer (without extention)	*/
#define	rwli	REC_WLIQ_GLACIER	/*	   recover file for WLIQ in glacier (without extention)	*/
#define	rwii	REC_WICE_GLACIER	/*	   recover file for WICE in glacier (without extention)	*/
#define	rTi	REC_GLACIER_TEMPERATURE	/*	   recover file for glacier temperature (without extention)	*/
#define	rns	REC_NUMBER_SNOWLAYER	/*	      recover file for number of snow layers(without extention)	*/
#define	rni	REC_NUMBER_GLACIERLAYER	/*	     recover file for number of glacier layers(without extention)	*/
#define	rsnag	REC_DIMENSIONLESS_SNOW_AGE	/*	  recover file for dimensionless snow age (without extention)	*/
#define	rhsup	REC__h_sup	/*	 recover file for h_sup (without extention)	*/
#define	rwt	REC__wt	/*	 recover file for wt (without extention)	*/
#define	rQch	REC__Qchannel	/*	 recover file for Qchannel (without extention)	*/
#define	rSFA	REC_SFA	/*	 recover file for SFR (without extention)	*/
#define	fHpatch	zSFA	/*	 recover file for zSFR (without extention)	*/


