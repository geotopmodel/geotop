rm(list=ls())

soiltn <- '/home/ecor/local/src/geotop_dev/tests/3D/panola_25pixel_nobed_hydrostatic/input/soil-classes/soil0001.txt' 
soilto <- '/home/ecor/local/src/geotop_dev/tests/3D/panola_25pixel_nobed_hydrostatic/input/soil-classes/soil0001.old.txt' 
soildf <- read.table(soilto,header=TRUE,sep=",")
cosz <-  cos(14.238626/180*pi)
dz <- soildf$Dz
z <- dz
for (i in 2:length(z)) {
	
	z[i] <- z[i-1]+dz[i]
	
}

z <- z-dz/2
d <- max(z)

soildf$InitPsi <- (z-d)*cosz
write.table(soildf,file=soiltn,row.names=FALSE,sep=",",quote=FALSE)