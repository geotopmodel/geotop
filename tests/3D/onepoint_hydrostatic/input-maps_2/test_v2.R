rm(list=ls())

library(geotopbricks)

from <- '/home/ecor/local/src/geotop_dev/tests/3D/onepoint_hydrostatic/input-maps_2' 
to <- '/home/ecor/local/src/geotop_dev/tests/3D/onepoint_hydrostatic/input-maps' 


maps <- list.files(from,pattern=".asc",full.name=TRUE)
names(maps) <- list.files(from,pattern=".asc",full.name=FALSE)
newmaps <- paste(to,names(maps),sep="/")
names(newmaps) <- names(maps)

rmap <- lapply(X=maps,FUN=raster)
rmapn <- rmap

for (it in names(maps)) {
	
	val <- rmapn[[it]][5,5]
	rmapn[[it]][] <- NA
	rmapn[[it]][5,3:5] <- val
	writeRasterxGEOtop(x=rmapn[[it]],filename=newmaps[[it]])
	
}