# TODO: Add comment
# 
# Author: stgruber
###############################################################################

#switch
task<-"PREP" #"PROC" #"PREP"

#source gt interaction
source("gt_control.R")


#==============================================================================
# MAKE EXPERIMENT STEP 1: establish sensitivity to soil and time discretization
# -	time stepping (2 versions, hourly and reduced)
# -	soil depth
# -	soil resolution
# - long spin-up
#==============================================================================
gt_example.01<-function(){
#paths	
eroot_loc<-"/Users/stgruber/Desktop/sonnblick/example01"
eroot_sim<-"/group/geotop/sim/permanet/sonnblick/example01"

#Schršdiger settings
eroot_sch<-"/lustre/stgruber/gt_new/example01"
name_sch<-"geotop_ex01"
maxtime_sch<-"3600"   #max time [s] per CPU

#time frame
InitDateDDMMYYYYhhmm   = "01/09/1988 00:00, 01/09/1988 00:00"
EndDateDDMMYYYYhhmm    = "31/08/1993 23:00, 31/08/1993 23:00"
NumSimulationTimes     = "20, 1"
DtPlotPoint            = "240, 3"

#sensitivity experiments
TimeStepEnergyAndWater = 3600

#read master directory
fs<-gt.exp.rmaster(eroot_loc)

#define soil  (rock,debris)
ThermalConductivitySoilSolids = 2.5
ThermalCapacitySoilSolids     = 1.E6

#define topography points and horizon files
hor_seq<-seq(from=   0,to=   0,by=  0)
ele_seq<-seq(from=2000,to=3250,by=250) # [m]
slp_seq<-seq(from=   0,to=  60,by= 20) # [deg]
asp_seq<-seq(from=   0,to= 270,by= 90) # [deg]
soi_seq<-seq(from=   1,to=   1,by=  1)
lan_seq<-seq(from=   1,to=   1,by=  1) 
dis_seq<-seq(from=  10,to=  10,by= 10) # [m] 
fre_seq<-seq(from=   1,to=   1,by=  1) # [mm]
swe_seq<-seq(from=1000000, to=1000000,  by=1) # [mm]
cur_seq<-seq(from=0, to=0,  by=0) # [1/m]

#make standard horizon files
horfile.write(eroot_loc)
#make points and get data frame
topo<-points.make(eroot_loc, hor_seq, ele_seq, slp_seq, asp_seq, lan_seq,soi_seq, dis_seq, fre_seq, swe_seq, cur_seq)

#set simulation period and timing			
fs<-gt.par.wvar(fs,"TimeStepEnergyAndWater",TimeStepEnergyAndWater[timing])	
fs<-gt.par.wvar(fs,"InitDateDDMMYYYYhhmm",InitDateDDMMYYYYhhmm)			
fs<-gt.par.wvar(fs,"EndDateDDMMYYYYhhmm" ,EndDateDDMMYYYYhhmm)
fs<-gt.par.wvar(fs,"NumSimulationTimes"  ,NumSimulationTimes)
fs<-gt.par.wvar(fs,"DtPlotPoint"         ,DtPlotPoint)

#other fixed stuff
fs<-gt.par.wvar(fs,"InitGlacierDepth",0)

#change snow amount
fs<-gt.par.wvar(fs,"SnowCorrFactor",1)
fs<-gt.par.wvar(fs,"RainCorrFactor",1)

out<-NULL
enumber<-1 #initial experiment numer
counter<-1 #counter, always start at 1
#loop over experiments ()
	for (soila in seq(from=1,to=3, by=1)) {
		for (soilb in seq(from=1,to=3, by=1)) {
			for (soilc in seq(from=1,to=3, by=1)) {
				for (soild in seq(from=1,to=2, by=1)) {		
					#get soil discretization
					soil.disc<-gt.sf.dz(dzmin[soila],zmax[soilb],base[soilc])
					#change soil discretization
					fs<-gt.par.wvar(fs,"SoilLayerThicknesses",soil.disc$dz)
					fs<-gt.par.wvar(fs,"SoilLayerNumber",length(soil.disc$dz)) 
					#make soil type
					fs<-gt.par.wvar(fs,"NormalHydrConductivity",NormalHydrConductivity[soild]) 
					fs<-gt.par.wvar(fs,"LateralHydrConductivity",LateralHydrConductivity[soild])
					fs<-gt.par.wvar(fs,"ThetaRes",ThetaRes[soild])               
					fs<-gt.par.wvar(fs,"ThetaSat",ThetaSat[soild])
					fs<-gt.par.wvar(fs,"AlphaVanGenuchten",AlphaVanGenuchten[soild])
					fs<-gt.par.wvar(fs,"NVanGenuchten",NVanGenuchten[soild])
					fs<-gt.par.wvar(fs,"VMualem",VMualem[soild])
					fs<-gt.par.wvar(fs,"SpecificStorativity",SpecificStorativity[soild])
					
					#switch purpose
					if (task == "PREP") {
						#write experiment
						epath<-gt.exp.write(eroot_loc,eroot_sim,enumber,fs)
						#write points and increment
						write.csv(topo, file = epath&"/listpoints.txt", quote=FALSE,row.names = FALSE)
					} else {
						list<-exp.read.lowest.annual(eroot_loc,enumber)
						list$ThataSat<-ThetaSat[soild]
						list$dzmin<-dzmin[soila]
						list$zmax<-zmax[soilb]
						list$base<-base[soilc]
						list$timing<-timing
						list$enumber<-enumber
					}	
					enumber<-enumber+1 #increment experiment number
					counter<-counter+1 #increment experiment number
					out<-rbind(out,list)
					write.table(out, file = eroot_loc&"/summary_exp1.csv")
					print(enumber)
}}}}
#write Schroedinger Batch file for SUN gridengine
gt.exp.schroedinger(eroot_loc, eroot_sch, name_sch, maxtime_sch)

}
#-- finished experiment one ------------------

