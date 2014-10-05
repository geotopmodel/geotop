# TODO: Add comment
# 
# Author: stgruber
###############################################################################



#==============================================================================
# MISC STUFF
#==============================================================================
#construct method to concatenate strings with ampersand
"&" <- function(...) UseMethod("&") 
"&.default" <- .Primitive("&") 
"&.character" <- function(...) paste(...,sep="") 

require("gregmisc")



#==============================================================================
# FIND LINE WITH A SPECIFIC KEYWORD (CASE INSENSITIVE)
#==============================================================================
gt.par.fline<-function(fs, keyword) { 
	nl<-length(fs)
	comchar<-"!" #character to indicate comments
	splchar<-"=" #sparate variable name and content
	keyword<-toupper(keyword) #make case insensitive
	found<-FALSE
	
	#loop
	for (l in 1:nl) {
		line<-trim(fs[l]) #assign line, but first remove leading/trailing blanks
		if ((substr(line,1,1) != comchar)*(nchar(line) >= 1)) {
			varname<-unlist(strsplit(line,splchar, fixed = TRUE))[1]
			varname<-trim(varname) #remove blanks
			if ((toupper(varname) == keyword)) found<-TRUE
		}
		if (found == TRUE) break #exit loop if keyword was found
	}	
	
	#return result
	if (found == FALSE) {
		return(NA)
	} else {
		return(l)
	}
}


#==============================================================================
# EXTRACT VARIABLE IN A SPECIFIC LINE
#==============================================================================
gt.par.rline<-function(fs,ln) { 
	#fs --> file string containing entire parameter file
	#ln --> line number to read
	comchar<-"!" #indicate comments
	splchar<-"=" #sparate variable name and content
	vecchar<-"," #separate vector elements

	line<-trim(fs[ln]) #assign line, but first remove leading/trailing blanks
	line<-unlist(strsplit(line,splchar, fixed = TRUE))
	varname<-line[1] #get variable name
	varvalu<-lnie[2] #get variable value
	varvalu<-unlist(strsplit(line,comchar, fixed = TRUE))[1] #delete comment at end of line if one exists	
	varvalu<-unlist(strsplit(line,vecchar, fixed = TRUE))    #separate vector elements
	varvalu<-trim(varvalu) #remove blanks
	
	#return result
	return(varvalu)
}


#==============================================================================
# WRITE VARIABLE INTO A SPECIFIC LINE
#==============================================================================
gt.par.wline<-function(fs,ln,vs) { 
	#fs --> file string containing entire parameter file
	#ln --> line number to read
	#vs --> variable character string to insert
	comchar<-"!" #indicate comments
	splchar<-"=" #sparate variable name and content
	vecchar<-"," #separate vector elements

	#check is character string is vector or scalar
	if (length(vs) > 1) vs<-paste(vs,sep="",collapse=vecchar)
	line<-trim(fs[ln]) #assign line, but first remove leading/trailing blanks
	varname<-unlist(strsplit(line,splchar,fixed = TRUE))[1] #get variable name	
	line<-varname&splchar&vs
	fs[ln]<-line
	return(fs)
}


#==============================================================================
# WRITE VARIABLE INTO FILESTRING OF PARAMETER FILE
#==============================================================================
gt.par.wvar<-function(fs,keyword,vs,type="NUMERIC") { 
	#fs --> file string containing entire parameter file
	#ln --> line number to read
	#vs --> variable to insert
	comchar<-"!" #indicate comments
	splchar<-"=" #sparate variable name and content
	vecchar<-"," #separate vector elements
	
	if (toupper(type) == "STRING") {
		vs<-"\""&vs&"\""
	} else {
		#convert to character
		vs<-as.character(vs)	
	}
	
	#check is character string is vector or scalar
	if (length(vs) > 1) vs<-paste(vs,sep="",collapse=vecchar)
	
	#find line
	ln<-gt.par.fline(fs, keyword)

	if (is.na(ln) == FALSE) { #insert into existing line
		fs<-gt.par.wline(fs,ln,vs)
	} else {        #if line does not exist, append at end
		line<-keyword&splchar&vs
		fs<-c(fs,line)
	}
	return(fs)
}


#==============================================================================
# READ EXPERIMENT MASTER PARAMETER FILE
#==============================================================================
#generate experiment directory and links
gt.exp.rmaster<-function(eroot) {
	#DIRECTORY STRUCTURE	
	#experiment_root/           --> base directory of an experiment
	#experiment_root/_master/   --> simulation template, experiments use common data from here
	#experiment_root/_control/  --> executable, source, batch.txt, documentation
	#experiment_root/000001/    --> path for experiment 1
	#fs 						--> filestring with modified parameter file to define experiment
	parfilename<-"geotop.inpts" #standard name of parameter file
	fs<-readLines(eroot&"/_master/"&parfilename) 
	return(fs)
}


#==============================================================================
# WRITE EXPERIMENT DIRECTORY AND PARAMETER FILE
#==============================================================================
#generate experiment directory and links
gt.exp.write<-function(eroot_loc,eroot_sim,enumber,fs) {
	#DIRECTORY STRUCTURE	
	#experiment_root/           --> base directory of an experiment
	#experiment_root/_master/   --> simulation template, experiments use common data from here
	#experiment_root/_control/  --> executable, source, batch.txt, documentation
	#experiment_root/000001/    --> path for experiment 1
	#fs 						--> filestring with modified parameter file to define experiment
	#eroot_loc                  --> experiment root on local machine (where simulation is prepared)
	#eroot_sim                  --> experiment root on simulation machine (used for batch file writing)
	comchar<-"!" #indicate comments
	parfilename<-"geotop.inpts" #standard name of parameter file
	#make directory ------------------------------------
	enumber<-round(enumber,0)
	epath<-eroot_loc&"/"&formatC(enumber,width=6,flag="0")
	dir.create(epath, showWarnings = TRUE, recursive = FALSE)

	#write paramater_file
	con <- file(epath&"/"&parfilename, "w")  # open an output file connection
	cat(comchar,"SCRIPT-GENERATED EXPERIMENT FILE",'\n', file = con,sep="")
	cat(fs, file = con,sep='\n')
	close(con)
	
	#append run call to batch file ---------------------     
	# put "#!/bin/sh"   or "#!/usr/bin/env tcsh" into first line
	batchfile_normal<-eroot_loc&"/_control/batch.txt"
	out_normal <- file(batchfile_normal, "a")
	epath_sim<-eroot_sim&"/"&formatC(enumber,width=6,flag="0")
	#cat(eroot_sim,"/_control/GEOtop ", epath_sim,"/&","\n",file=out,sep="")
	cat(eroot_sim,"/_control/GEOtop ", epath_sim,"\n",file=out_normal,sep="")
	close(out_normal)
	return(epath)
}

#==============================================================================
# WRITE SCHROEDINGER BATCH FILE
#==============================================================================
# this script will make a list of experiment directories contained in a
# local simmulation directory (eroot_loc) and write a bach file for the
# Schroedinger gridengine. This will assume the file geotop_in.tar.gz to be
# present in the Schroedinger simmulation directory (eroot_loc) that contains
# all input directories. This file will be unpacked, simulations performed and
# afterwards, simulation results will be compressed into geotop_out.tar.gz.
# Execute on Schroedinger by typing "qsub myscript.sh" (or any other filename)
# You can control the state of your job by "qstat -u $USER"
gt.exp.schroedinger<-function(eroot_loc, eroot_sim, eroot_sch, name, maxtime) {
	#settings
	batchfilename<-"batch.txt" #standard name of parameter file
	batchfilename_schroedinger <-"batch_schroedinger.txt" #file with list of individual simulations
	scriptfilename_schroedinger<-"geotop_schroedinger.sh" #script for job submission: "qsub myscript.sh" (or any other filename)

	#---- make Schroedinger BATCH file ----------------------------------------
	#get list of 
	fs<-readLines(eroot_loc&"/_control/"&batchfilename) 

	#count lines, change directories
	totalnumber<-length(fs) #number of lines
	fs<-gsub(eroot_sim, eroot_sch, fs,ignore.case = FALSE, fixed = TRUE)
	fs<-gsub("&", "", fs,ignore.case = FALSE, fixed = TRUE)
	
	#write SChroedinger batch file
	con <- file(eroot_loc&"/_control/"&batchfilename_schroedinger, "w")  # open an output file connection
	cat(fs,file = con,sep='\n')
	close(con)
	
	
	#---- make Schroedinger SCRIPT file ---------------------------------------
	#Create filestring and insert first line, Bourne shell
	fs<-"#!/bin/sh"	
	
	#set number of simulations and step size
	fs<-c(fs,"#$ -t 1-"&round(totalnumber,0)&":8      # setting range and step size")
	
	#set maximum time in seconds	
	fs<-c(fs,"#$ -N="&as.character(name)&" # name of this Schroedinger run")
	
	#set maximum time in seconds	
	fs<-c(fs,"#$ -l s_cpu="&round(maxtime,0)&" # max time [s] per cpu")
	
    #diverse settings
	fs<-c(fs,"#$ -S /bin/sh     # shell used")
	fs<-c(fs,"#$ -v PATH        # environment variables to be exported to the execution context of the job")
	fs<-c(fs,"#$ -o $JOB_NAME_$JOB_ID_$TASK_ID.out  # create output file per task ($TASK_ID only valid here)")
	fs<-c(fs,"#$ -j y           # error stream of job merged into standard output stream")
	fs<-c(fs,"# this script launches 8 simulations on one node and the startes the next 8")
	
	# set main variables
	fs<-c(fs,"eroot_sch="&eroot_sch)
	
	#entering loop 0 to 7
	fs<-c(fs,"#entering loop 0 to 7")	
	fs<-c(fs,"for fi in 0 1 2 3 4 5 6 7")	
	fs<-c(fs,"do")	
	fs<-c(fs,"  counter=$(( $SGE_TASK_ID + $fi ))   # increase counter, bash arithmetics, here SGE_TASK_ID not as above")	
	#use paste to insert quotations into string
	fs<-c(fs,paste("  gtjob=`sed -n ","\"","${counter}p","\""," $eroot_sch/_control/batch_schroedinger.txt` #get line from job list",sep=""))
	fs<-c(fs,"  echo 'counter=$counter, $gtjob' # provide terminal feedback of what is done") 
	fs<-c(fs,"  $gtjob & # execute job in background --> comment this line for testing")
	fs<-c(fs,"done")
	
	#wait for jobs to be done
	fs<-c(fs,"wait # wait until all background jobs are finished")
	
	#write file
	con <- file(eroot_loc&"/_control/"&scriptfilename_schroedinger, "w")  # open an output file connection
	cat(fs,file = con,sep='\n')
	close(con)
}



#==============================================================================
# MAKE HORIZON FILES
#==============================================================================
#Needs to be updated
#make series of standard horizon files - e.g.,  hor0001.txt for 10¡ 
horfile.write<-function(eroot_loc) {
	for (el in seq(from=0,to=90,by=10)) {
		hor<-data.frame(az=c(0,360),el=c(el,el))
		#write file
		outfile<-eroot_loc&"/_master/hor000"&round(el/10+1,0)&".txt"
		write.csv(hor,file=outfile,quote=FALSE,row.names = FALSE)
	}
}
#==============================================================================


#==============================================================================
# SKY VIEW FACTOR (analtic for slope and horizon of constant elevation, only)
#==============================================================================
sky.view<-function(slp,hor.el) {
	sky<-cos((slp+hor.el)/2*pi/180)^2
	return(round(sky,2))
}
#==============================================================================


#==============================================================================
# MAKE POINTS FOR POINT FILES
#==============================================================================
points.make<-function(eroot_loc, hor_seq, ele_seq, slp_seq, asp_seq, lan_seq, 
		              soi_seq, dis_seq, fre_seq, swe_seq, cur_seq) {	
	#initialize names
	ID_vec<-NULL
ele_vec<-NULL
slp_vec<-NULL
asp_vec<-NULL
sky_vec<-NULL
hor_vec<-NULL #horizon elevation (uniform)	
lan_vec<-NULL #land cover type number
soi_vec<-NULL #soil type number
dis_vec<-NULL #distance to free drainage
fre_vec<-NULL #depression of free drainage below surface
swe_vec<-NULL #maximum snow water equivalent
cur_vec<-NULL #curvature (snow transport) in all directions

#make loop and generate data frame
ID<-1
for (ele in ele_seq) {
	for (slp in slp_seq) {
		if (slp == 0) {avec<-0:0} else {avec<-asp_seq}
		for (asp in avec) {
			for (hor in hor_seq) {
				sky<-sky.view(slp,hor) 
				for (lan in lan_seq) {					
					for (soi in soi_seq) {
						for (dis in dis_seq) {
							for (fre in fre_seq) {
								for (swe in swe_seq) {
									for (cur in cur_seq) {
										ID_vec<-c( ID_vec,ID)
										ele_vec<-c(ele_vec,ele)
										slp_vec<-c(slp_vec,slp)
										asp_vec<-c(asp_vec,asp)
										sky_vec<-c(sky_vec,sky)
										hor_vec<-c(hor_vec,round(hor/10+1,0)) #make right number		
										lan_vec<-c(lan_vec,lan)
										soi_vec<-c(soi_vec,soi)
										dis_vec<-c(dis_vec,dis)
										fre_vec<-c(fre_vec,fre)
										tmp<-swe*(90-slp)
										if (slp <= 50) {tmp<-1000000}
										swe_vec<-c(swe_vec,tmp)
										cur_vec<-c(cur_vec,cur)		
										ID<-ID+1
									}
								}
							}
						}
					}
				}	
			}	
		}	
	}	
}	
	topo<-data.frame(ID=ID_vec,ele=ele_vec,slp=slp_vec,asp=asp_vec,sky=sky_vec,
			landcover=lan_vec,soil=soi_vec,dist=dis_vec,free=fre_vec, 
			maxswe=swe_vec,curvNS=cur_vec,curvWE=cur_vec,
			curvNwSe=cur_vec,curvNeSw=cur_vec,hor=hor_vec)
	return(topo)
}
#==============================================================================



#==============================================================================
#FUNCTION TO MAKE THICKNESS AND DEPTH OF LAYERS
#============================================================================== 
#  dzmin: minimal z-spacing ("how fine is the grid") 
#   zmax: depth of lowermost node center ("how large is the domain") 
#   base: resolution reduction ("how strong is the grid corsened with depth")
gt.sf.dz <- function(dzmin,zmax,base) {
	#layer thicknesses	
	dz<-dzmin*base^(0:2005)
	
	#depth of layer lower boundary
	z_bot<-cumsum(dz)
	#depth of layer upper boundary
	z_top<-z_bot-dz
	#depth of layer center (node)
	z<-(z_bot+z_top)/2
	
	#data frame
	discr<-data.frame(dz=dz,z=z,z_bot=z_bot,z_top=z_top)
	
	#restrict to maximum depth
	discr<-subset(discr,z < zmax)
	nz<-length(discr$z)
	
	discr$dz[nz] <- zmax-discr$z_bot[nz-1] 	
	discr$z[nz]  <- discr$z_bot[nz-1] + discr$dz[nz]/2    
	
	#drop auxiliary columns
	discr<-discr[1:nz,1:2] 
    return(discr)
}


#==============================================================================
#FUNCTION TO MAKE SNOW PARAMETERIZATIONS
#============================================================================== 
#  dzmin: minimal z-spacing ("how fine is the grid") 
#   zmax: depth of lowermost node center ("how large is the domain") 
#   base: resolution reduction ("how strong is the grid corsened with depth")
#         (has the same mechanism as in soil)
#   nlow: number of finely discretized layers near the soil 
#         (the layer above has near-infinite maximum thickness)
#   The maximum thickness of a layer is the sum of its minimum thickness and that
#   of the next thicker layer.
gt.snow.dz <- function(dzmin,zmax,base,nlow) {
	#bottom layer thicknesses	
	dzmin_bot<-dzmin*base^(0:(nlow-1))
	dzmax_bot<-dzmin_bot+dzmin*base^(1:nlow) #max is equal to sum of min of layer and that above
	
    #check dimensions. stop function if lower nodes take too much room
	z_bot<-max(cumsum(dzmax_bot))
	if (zmax/z_bot <= 2) {return(NA)}
	
	#top layer thicknesses	
	dzmin_top<-dzmin*base^(500:0)
	dzmax_top<-dzmin_top+dzmin*base^(501:1) #max is equal to sum of min of layer and that above
	
	#restrict to maximum thickness
	ind<-rev(cumsum(rev(dzmax_top)))<=(zmax-z_bot)
	
	#combine
	dzmi<-c(dzmin_bot,dzmin_top[ind])
	dzma<-c(dzmax_bot,dzmax_top[ind])
	
	#assign infinite layer
	dzma[nlow+1]<-1000000
	#return data frame
	return(data.frame(dzmin=round(dzmi,2),dzmax=round(dzma,2)))
}



#==============================================================================
# ACCESS EXPERIMENT RESULTS
#==============================================================================
read.run<-function(eroot,enumber,file){
	#make directory
	path<-eroot&"/"&formatC(enumber,width=6,flag="0")&"/"
	#read data
	return(read.csv(path&file))
}




#==============================================================================
# READ GEOTOP MATRIX INTO DATA FRAME 
#============================================================================== 
#usage to extract comment:
#comment<-gt.rmatrix("/group/geotop/sometestfile.txt",comment=TRUE)
#usage to read data:
#comment<-gt.rmatrix("/group/geotop/sometestfile.txt")
gt.rmatrix<-function(infile,comment=FALSE) {
	comchar<-"!" #character to indicate comments
	GTtime<-"%d/%m/%Y %H:%M" #descrition of data format used by GEOtop
	#count comment lines
	found<-FALSE #begin of data found? (first line without comment)
	fs <-readLines(infile,n=1000)
	for (l in 1:length(fs)) {
		line<-trim(fs[l])
		if (substr(line,1,1) != comchar) break
	}
	l<-l-1
	#decide if comment of data is returned
	if (comment == TRUE) {
		return(fs[1:l]) #get and return comment
	} else {
		data<-read.csv(infile, skip=l, header = TRUE, sep = ",", dec=".")
		#if first column is a date, then convert
		#data format always looks like this:
		#31/01/2004 13:24:00
		#1234567890123456789
		isDate<-TRUE
		data[,1]<-as.character(data[,1])
		if (nchar(data[1,1]) != 16)        isDate<-FALSE
		if (substr(data[1,1], 3, 3) != "/") isDate<-FALSE
		if (substr(data[1,1], 6, 6) != "/") isDate<-FALSE
		if (substr(data[1,1],11,11) != " ") isDate<-FALSE
		if (substr(data[1,1],14,14) != ":") isDate<-FALSE	
		if (isDate == TRUE) data$Date <- as.POSIXct(data[,1], GTtime, tz="UTC")	
		return(data)
	}
}		



#==============================================================================
# WRITE DATA FRAME INTO GEOTOP MATRIX  
#============================================================================== 
gt.wmatrix<-function(data, outfile, comment="") {
	comchar<-"!" #character to indicate comments
	GTtime<-"%d/%m/%Y %H:%M:%S" #descrition of data format used by GEOtop
	#count comment lines
	found<-FALSE #begin of data found? (first line without comment)
	fs <-readLines(infile,n=1000)
	for (l in 1:length(fs)) {
		line<-trim(fs[l])
		if (substr(line,1,1) != comchar) break
	}
	comment<-fs[1:l] #get comment
	data<-read.csv(infile, skip=l-1, header = TRUE, sep = ",", dec=".")
	
	#if first column is a date, then convert
	#data format always looks like this:
	#31/01/2004 13:24:00
	#1234567890123456789
	isDate<-TRUE
	if (nchar(data[1,1]) != 19)        isDate<-FALSE
	if (substr(data[1,1],3, 1) != "/") isDate<-FALSE
	if (substr(data[1,1],6, 1) != "/") isDate<-FALSE	
	if (substr(data[1,1],11,1) != " ") isDate<-FALSE
	if (substr(data[1,1],14,1) != ":") isDate<-FALSE	
	if (substr(data[1,1],17,1) != ":") isDate<-FALSE
	if (isDate == TRUE) data$Date <- as.POSIXct("31/01/2004 13:24:00", GTtime, tz="UTC")
	res<-list(comment=comment, data=data)	
	return(res)
}	


