#include "stdio.h"
#include "math.h"

void main() {

FILE *f;
int nc=11, nr=9, c, r, y, m, d, h, days_in_month[12];
double ds=500.0, jd=734868.999999, prec=10.0, conc=0.001, temp=30.0;
double xll=600000.0, yll=5000000.0, no_data=-9999.0;


days_in_month[0]=31; // Januar
days_in_month[1]=28; // Februar
days_in_month[2]=31; // Maerz
days_in_month[3]=30; // April
days_in_month[4]=31; // Mai
days_in_month[5]=30; // Juni
days_in_month[6]=31; // Juli
days_in_month[7]=31; // August
days_in_month[8]=30; // September
days_in_month[9]=31; // Oktober
days_in_month[10]=30; // November
days_in_month[11]=31; // Dezember


f=fopen("dem.asc","w");
fprintf(f,"ncols %i\nnrows %i\nxllcorner %f\nyllcorner %f\ncellsize %f\nNODATA_value %f\n",nc,nr,xll,yll,ds,no_data);
for (r=1; r<=nr; r++){
	for (c=1; c<=nc; c++){
		fprintf(f,"%f ",2001.0 - (double)c );
	}
	fprintf(f,"\n");
}
fclose(f);

f=fopen("net.asc","w");
fprintf(f,"ncols\t\t%i\nnrows\t\t%i\nxllcorner\t600000.0\nyllcorner\t5000000.0\ncellsize\t%f\nNODATA_value\t-9999.0\n",nc,nr,ds);
for (r=1; r<=nr; r++){
	for (c=1; c<=nc; c++){
		if (r==5 && c>3){fprintf(f,"%f\t",10.0);}
		else {fprintf(f,"%f\t",0.0);}
	}
	fprintf(f,"\n");
}
fclose(f);

/*
f=fopen("surface.asc","w");
fprintf(f,"ncols\t\t%i\nnrows\t\t%i\nxllcorner\t600000.0\nyllcorner\t5000000.0\ncellsize\t%f\nNODATA_value\t-9999.0\n",nc,nr,ds);
for (r=1; r<=nr; r++){
	for (c=1; c<=nc; c++){
		if (r==5 && c>3){fprintf(f,"%f\t",40.0);}
		else {fprintf(f,"%f\t",20.0);}
	}
	fprintf(f,"\n");
}
fclose(f);

f=fopen("landcover.asc","w");
fprintf(f,"ncols\t\t%i\nnrows\t\t%i\nxllcorner\t600000.0\nyllcorner\t5000000.0\ncellsize\t%f\nNODATA_value\t-9999.0\n",nc,nr,ds);
for (r=1; r<=nr; r++){
	for (c=1; c<=nc; c++){
		fprintf(f,"%f\t",6.0);
		
	}
	fprintf(f,"\n");
}
fclose(f);

f=fopen("watertable.asc","w");
fprintf(f,"ncols\t\t%i\nnrows\t\t%i\nxllcorner\t600000.0\nyllcorner\t5000000.0\ncellsize\t%f\nNODATA_value\t-9999.0\n",nc,nr,ds);
for (r=1; r<=nr; r++){
	for (c=1; c<=nc; c++){
		fprintf(f,"%f\t",-1999.0 - (double)c );
	}
	fprintf(f,"\n");
}
fclose(f);
*/

f=fopen("meteo0001.txt","w");
fprintf(f,"Date,JDfrom0,Iprec,AirT,Conc\n");
for (y=2012; y<=2013; y++){
	for (m=1; m<=12; m++){
		for (d=1; d<=days_in_month[m-1]; d++){
			for (h=0; h<24; h++){
				if (m<10){
					if (d<10){
						if (h<10) {fprintf(f,"0%i/0%i/%i 0%i:00,%f,",d,m,y,h,jd);}
						else {fprintf(f,"0%i/0%i/%i %i:00,%f,",d,m,y,h,jd);}
					}
					else {
						if (h<10) {fprintf(f,"%i/0%i/%i 0%i:00,%f,",d,m,y,h,jd);}
						else {fprintf(f,"%i/0%i/%i %i:00,%f,",d,m,y,h,jd);}
					}
				}
				else {
					if (d<10){
						if (h<10) {fprintf(f,"0%i/%i/%i 0%i:00,%f,",d,m,y,h,jd);}
						else {fprintf(f,"0%i/%i/%i %i:00,%f,",d,m,y,h,jd);}
					}
					else {
						if (h<10) {fprintf(f,"%i/%i/%i 0%i:00,%f,",d,m,y,h,jd);}
						else {fprintf(f,"%i/%i/%i %i:00,%f,",d,m,y,h,jd);}
					}
				}
				fprintf(f,"%f,%f,%f\n",prec,temp,conc);
				if (jd>735384.999999+ 20.0-2.0/24.0){prec=0.0; temp=10.0;}
				if (jd>735384.999999+ 40.0-2.0/24.0){prec=5.0; temp=30.0;}
				if (jd>735384.999999+ 60.0-2.0/24.0){prec=0.0; temp=10.0;}
				if (jd>735384.999999+ 80.0-2.0/24.0){prec=2.5; temp=30.0;}
				if (jd>735384.999999+100.0-2.0/24.0){prec=0.0; temp=10.0;}
				jd=jd+1.0/24.0;
			}
		}
	}
}
fclose(f);


f=fopen("geotop.inpts","w");
fprintf(f,"!=======================================\n");
fprintf(f,"! INPUT FOR GEOTOP V2.0\n");
fprintf(f,"! atrifical input\n");
fprintf(f,"!=======================================\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! GENERAL SETTINGS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"! Calculation max time step dt in s\n");
fprintf(f,"TimeStepEnergyAndWater\t\t\t=\t3600\n");
fprintf(f,"InitDateDDMMYYYYhhmm\t\t\t=\t01/05/2013 00:00\n");
fprintf(f,"EndDateDDMMYYYYhhmm\t\t\t=\t10/05/2013 00:00\n\n");

fprintf(f,"! Catchment centroid (for Sun position)\n");
fprintf(f,"Latitude\t\t\t\t=\t46.75\n");
fprintf(f,"Longitude\t\t\t\t=\t10.70\n");
fprintf(f,"StandardTimeSimulation\t\t\t=\t1\n\n");

fprintf(f,"! Simulation settings\n");
fprintf(f,"WaterBalance\t\t\t\t=\t1\n");
fprintf(f,"EnergyBalance\t\t\t\t=\t1\n");
fprintf(f,"TransportModel\t\t\t\t=\t1\n\n");

fprintf(f,"! Output timeseries Dt in hours\n");
fprintf(f,"DtPlotDischarge\t\t\t\t=\t1\n");
fprintf(f,"DtPlotPoint\t\t\t\t=\t1\n");
fprintf(f,"DtPlotBasin\t\t\t\t=\t24\n\n");

fprintf(f,"! Output maps Dt in hours\n");
fprintf(f,"OutputSoilMaps\t\t\t\t=\t24\n");
fprintf(f,"OutputSurfEBALMaps\t\t\t=\t24\n");
fprintf(f,"OutputMeteoMaps\t\t\t\t=\t24\n");
fprintf(f,"OutputSnowMaps\t\t\t\t=\t24\n");
fprintf(f,"OutputGlacierMaps\t\t\t=\t24\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! METEO STATIONS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"NumberOfMeteoStations\t\t\t=\t1\n");
fprintf(f,"MeteoStationCoordinateX\t\t\t=\t601000\n");
fprintf(f,"MeteoStationCoordinateY\t\t\t=\t5001000\n");
fprintf(f,"MeteoStationElevation\t\t\t=\t1998\n");
fprintf(f,"MeteoStationWindVelocitySensorHeight\t=\t10\n");
fprintf(f,"MeteoStationTemperatureSensorHeight\t=\t10\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! METEO HEADERS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"HeaderDateDDMMYYYYhhmmMeteo\t\t=\t\"Date\"\n");
fprintf(f,"HeaderJulianDayfrom0Meteo\t\t=\t\"JDfrom0\"\n");
fprintf(f,"HeaderIPrec\t\t\t\t=\t\"Iprec\"\n");
fprintf(f,"HeaderConcentration\t\t\t=\t\"Conc\"\n");
fprintf(f,"HeaderAirTemp\t\t\t\t=\t\"AirT\"\n\n\n");
/*
fprintf(f,"HeaderWindVelocity\t\t\t=\t\"WindSp\"\n");
fprintf(f,"HeaderWindDirection\t\t\t=\t\"WindDir\"\n");
fprintf(f,"HeaderWindX\t\t\t\t=\t\"WindX\"\n");
fprintf(f,"HeaderWindY\t\t\t\t=\t\"WindY\"\n");
fprintf(f,"HeaderRH\t\t\t\t=\t\"RelHum\"\n");
fprintf(f,"HeaderAirTemp\t\t\t\t=\t\"AirT\"\n");
fprintf(f,"HeaderAirPress\t\t\t\t=\t\"AirP\"\n");
fprintf(f,"HeaderSWglobal\t\t\t\t=\t\"SWglobal\"\n");
fprintf(f,"HeaderCloudSWTransmissivity\t\t=\t\"CloudTrans\"\n\n\n");
*/

fprintf(f,"!=======================================\n");
fprintf(f,"! OUTPUT POINT SETTINGS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"CoordinatePointX\t\t\t=\t604500\n");
fprintf(f,"CoordinatePointY\t\t\t=\t5002000\n\n\n");


fprintf(f,"!=======================================\n"); 
fprintf(f,"! RATES OF DECREASE WITH ELEVATION\n");
fprintf(f,"!=======================================\n"); 
fprintf(f,"! K/1000m\n");
fprintf(f,"LapseRateTemp\t\t\t\t=\t9.8\n");
fprintf(f,"LapseRateDewTemp\t\t\t=\t2.5\n");
fprintf(f,"! mm/1000m\n");
fprintf(f,"LapseRatePrec\t\t\t\t=\t-0.3\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! LAND COVER\n");
fprintf(f,"!=======================================\n");
fprintf(f,"! Types: 1 urban, 2 agriculture, 3 forest, 4 grassland, 5 pasture, 6 bare rocks, 7 bare soil, 8 glacier, 9 lake/marshland\n\n");

fprintf(f,"NumLandCoverTypes\t\t\t=\t9\n");
fprintf(f,"SoilRoughness\t\t\t\t=\t100\n");
fprintf(f,"ThresSnowSoilRough\t\t\t=\t100\n");
fprintf(f,"VegHeight\t\t\t\t=\t0,0,8000,200,0,0,0,0,0\n");
fprintf(f,"ThresSnowVegUp\t\t\t\t=\t50\n");
fprintf(f,"ThresSnowVegDown\t\t\t=\t10\n");
fprintf(f,"LSAI\t\t\t\t\t=\t0,0,4,2,0,0,0,0,0\n");
fprintf(f,"CanopyFraction\t\t\t\t=\t0,0,1,1,0,0,0,0,0\n");
fprintf(f,"DecayCoeffCanopy\t\t\t=\t2.5\n");
fprintf(f,"VegSnowBurying\t\t\t\t=\t1\n");
fprintf(f,"RootDepth\t\t\t\t=\t0,0,500,200,0,0,0,0,0\n");
fprintf(f,"MinStomatalRes\t\t\t\t=\t60\n");
fprintf(f,"VegReflectVis\t\t\t\t=\t0.1\n");
fprintf(f,"VegReflNIR\t\t\t\t=\t0.58\n");
fprintf(f,"VegTransVis\t\t\t\t=\t0.05\n");
fprintf(f,"VegTransNIR\t\t\t\t=\t0.25\n");
fprintf(f,"LeafAngles\t\t\t\t=\t0\n");
fprintf(f,"CanDensSurface\t\t\t\t=\t0.5\n");
fprintf(f,"SoilAlbVisDry\t\t\t\t=\t0.15\n");
fprintf(f,"SoilAlbNIRDry\t\t\t\t=\t0.25\n");
fprintf(f,"SoilAlbVisWet\t\t\t\t=\t0.15\n");
fprintf(f,"SoilAlbNIRWet\t\t\t\t=\t0.25\n");
fprintf(f,"SoilEmissiv\t\t\t\t=\t0.96\n");
fprintf(f,"SurFlowResLand\t\t\t\t=\t0.5\n");
fprintf(f,"SurFlowResExp\t\t\t\t=\t0.667\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! SOIL\n");
fprintf(f,"!=======================================\n");
fprintf(f,"! [mm] \n");
fprintf(f,"InitWaterTableDepth\t\t\t=\t2000\n");
fprintf(f,"SoilLayerThicknesses\t\t\t=\t30,100,200,500,1200,3000,5000,12000,25000\n");
fprintf(f,"! [C]\n");
fprintf(f,"InitSoilTemp\t\t\t\t=\t0.5\n");
fprintf(f,"! [mm/s]\n");
fprintf(f,"NormalHydrConductivity\t\t\t=\t0.1,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.0001\n");
fprintf(f,"LateralHydrConductivity\t\t\t=\t0.0025, 0.0025, 0.00125, 0.0005, 0.00025, 0.000125, 0.00005, 0.000025,0.000025\n");
fprintf(f,"! [-]\n");
fprintf(f,"ThetaRes\t\t\t\t=\t0.08\n");
fprintf(f,"WiltingPoint\t\t\t\t=\t0.15\n");
fprintf(f,"FieldCapacity\t\t\t\t=\t0.3\n");
fprintf(f,"ThetaSat\t\t\t\t=\t0.5\n");
fprintf(f,"NVanGenuchten\t\t\t\t=\t1.3\n");
fprintf(f,"! [mm^-1]\n");
fprintf(f,"AlphaVanGenuchten\t\t\t=\t0.004\n");
fprintf(f,"!\n");
fprintf(f,"ThermalConductivitySoilSolids\t\t=\t2.5\n");
fprintf(f,"!\n");
fprintf(f,"ThermalCapacitySoilSolids\t\t=\t2.30E+06\n");
fprintf(f,"! [mm^-1]\n");
fprintf(f,"SpecificStorativity\t\t\t=\t1.00E-05\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! CHANNELS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"RatioChannelWidthPixelWidth\t\t=\t0.1\n");
fprintf(f,"! in [mm]\n");
fprintf(f,"ChannelDepression\t\t\t=\t2000\n");
fprintf(f,"ThresWaterDepthLandDown\t\t\t=\t5\n");
fprintf(f,"ThresWaterDepthLandUp\t\t\t=\t50\n");
fprintf(f,"SurFlowResChannel\t\t\t=\t20\n");
fprintf(f,"ThresWaterDepthChannelUp\t\t=\t50\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! SNOW AND GLACIERS\n");
fprintf(f,"!=======================================\n");
fprintf(f,"NumMaxSnowLayers\t\t\t=\t5\n");
fprintf(f,"InfiniteSnowLayer\t\t\t=\t2\n"); 
fprintf(f,"MinLayerThicknessSnow\t\t\t=\t5,250,50,5,5\n");
fprintf(f,"MaxLayerThicknessSnow\t\t\t=\t50,1.10,500,100,10\n");
fprintf(f,"FreshSnowReflVis\t\t\t=\t0.91\n");
fprintf(f,"FreshSnowReflNIR\t\t\t=\t0.68\n");		
fprintf(f,"InitGlacierDensity\t\t\t=\t700\n");
fprintf(f,"InitGlacierTemp\t\t\t\t=\t-5\n");
fprintf(f,"NumMaxGlacierLayers\t\t\t=\t3\n");
fprintf(f,"MinLayerThicknessGlacier\t\t=\t4000,100,10\n");
fprintf(f,"MaxLayerThicknessGlacier\t\t=\t1.10,8000,100\n\n\n");


fprintf(f,"!=======================================\n");
fprintf(f,"! Numerical parameters\n");
fprintf(f,"!=======================================\n");
fprintf(f,"RichardTol\t\t\t\t=\t1.E-7\n");
fprintf(f,"MinLambdaWater\t\t\t\t=\t1.E-7\n");
fprintf(f,"RichardMaxIter\t\t\t\t=\t120\n");
fprintf(f,"MaxTimesHalvingTimeStepWater\t\t=\t20\n");
fprintf(f,"MaxCourantSupFlowLand\t\t\t=\t0.1\n");
fprintf(f,"MaxCourantSupFlowChannel\t\t=\t0.1\n");
fprintf(f,"MinSupWaterDepthLand\t\t\t=\t1\n");
fprintf(f,"MinDiffSupWaterDepthLandChannel\t\t=\t50\n");
fprintf(f,"MinTimeStepSupFlow\t\t\t=\t1\n");
fprintf(f,"HeatEqTol\t\t\t\t=\t1.E-4\n");
fprintf(f,"HeatEqMaxIter\t\t\t\t=\t200\n");
fprintf(f,"MaxTimesHalvingTimeStepEnergy\t\t=\t5\n");
fprintf(f,"CanopyMaxIter\t\t\t\t=\t3\n"); 
fprintf(f,"BusingerMaxIter\t\t\t\t=\t3\n"); 
fprintf(f,"TsMaxIter\t\t\t\t=\t3\n");
fprintf(f,"LocMaxIter\t\t\t\t=\t3\n\n"); 

fprintf(f,"! Morphological parameters\n");
fprintf(f,"SlopeWeight\t\t\t\t=\t1\n"); 
fprintf(f,"CurvatureWeight\t\t\t\t=\t100\n");
fprintf(f,"NumLowPassFilterOnDemForAll\t\t=\t2\n");
fprintf(f,"NumLowPassFilterOnDemForCurv\t\t=\t20\n\n");

fprintf(f,"! Energy budget settings\n");
fprintf(f,"FlagSkyViewFactor\t\t\t=\t1\n");
fprintf(f,"LWinParameterization\t\t\t=\t7\n"); 
fprintf(f,"MoninObukhov\t\t\t\t=\t2\n"); 
fprintf(f,"CanopyStabCorrection\t\t\t=\t1\n\n\n"); 


fprintf(f,"!=======================================\n");
fprintf(f,"!  FILE NAMES\n");
fprintf(f,"!=======================================\n\n");

fprintf(f,"! Input files\n\n");

fprintf(f,"DemFile\t\t\t\t\t=\t\"dem\"\n");
fprintf(f,"MeteoFile\t\t\t\t=\t\"meteo\"\n");
fprintf(f,"RiverNetwork\t\t\t\t=\t\"net\"\n\n");


fprintf(f,"! Output files\n\n");

fprintf(f,"! Tabs\n");
fprintf(f,"DischargeFile\t\t\t\t=\t\"tabs/discharge\"\n");
fprintf(f,"ConcentrationFile\t\t\t=\t\"tabs/concentration\"\n\n");

fprintf(f,"! Maps\n");
//fprintf(f,"WaterTableDepthMapFile\t\t\t=\t\"maps/watertable\"\n\n");


fclose(f);

}

