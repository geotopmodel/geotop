#creo la mask per le stazioni

library(raster)
library(rgdal)

#creo il raster
angoli<-readOGR(dsn = "/home/francesco/git/mountaineering/geotop/tests/PSQL_test/",layer ="mask")
angoli_raster<-raster(angoli)
#setto la risoluzione del raster
res(angoli_raster)<-1000
#assegno alla mappa raster tutti valori=1
for(i in 1:ncell(angoli_raster)){
      angoli_raster[i]<-1    
}  


setwd("/home/francesco/git/mountaineering/geotop/tests/PSQL_test/")
writeRaster(angoli_raster, "stazioni.asc", overwrite=TRUE)

#ritaglio la mappa dei parametri di interesse nel file raster to_be_cut

#usare list
to_be_cut<-raster()
cut<-crop(to_be_cut, extent(trim(angoli_raster)))
writeRaster(cut, "cut.asc", overwrite=TRUE)

