#creo la mask per le stazioni

library(raster)
library(rgdal)

# Devi dargli il percorso completo della mappa compreso il nome della mappa e l estensione (asc)
# gli devi dare il raster e il numero di cifre decimali
# Ricorda che al raster devi assegnare come nodatavalue -9999


#**************************************************************************************************
#****************************  FUN WRITE RASTER GEOTOP FRIENDLY *********************************
#**************************************************************************************************
writeRasterMNT <- function(pathMap,raster,dp){
#pathMap percorso della mappa raster che dovrà essere scritta
#raster oggetto Raster*
#numero di cifre decimali relative alla mappa
require(raster)
writeRaster(round(raster,dp),pathMap, format="ascii",overwrite=TRUE)
funChangeHeader(pathMap)
}


#**********************************************************************************************************************************
#****************************  FUN CHANGE HEADER FROM UPPER CASE LETTERS TO LOWER CASE LETTERS    *********************************
#**********************************************************************************************************************************
funChangeHeader <- function(pathMap){
#pathMap percorso della mappa raster che dovrà essere scritta
fs <- readLines(pathMap,warn=FALSE)
key_low <- tolower(fs[1:5])
fs[1:5]<-paste(key_low[1:5])
writeLines(fs,pathMap)
}



#carico la mask



mask_stazioni<-raster("mask_stazioni.asc")
asp<-raster("asp.asc")
aspect<-raster("aspect.asc")
slp<-raster("slp.asc")
lc<-raster("lc.asc")
st<-raster("st.asc")
wrf_sky<-raster("wrf_sky.asc")
svf<-raster("svf.asc")
dem<-raster("dem.asc")



asp_cut<-crop(asp,extent(trim(mask_stazioni)))
aspect_cut<-crop(aspect,extent(trim(mask_stazioni)))
slp_cut<-crop(slp,extent(trim(mask_stazioni)))
lc_cut<-crop(lc,extent(trim(mask_stazioni)))
st_cut<-crop(st,extent(trim(mask_stazioni)))
wrf_sky_cut<-crop(wrf_sky,extent(trim(mask_stazioni)))
svf_cut<-crop(svf,extent(trim(mask_stazioni)))
dem_cut<-crop(dem,extent(trim(mask_stazioni)))


#set NA values to -9999
fun <- function(x) { x[is.na(x)] <- -9999; return(x)}
asp_cut2 <- calc(asp_cut, fun)
aspect_cut2 <- calc(aspect_cut, fun)
slp_cut2 <- calc(slp_cut, fun)
lc_cut2 <- calc(lc_cut, fun)
st_cut2 <- calc(st_cut, fun)
wrf_sky_cut2 <- calc(wrf_sky_cut, fun)
svf_cut2 <- calc(svf_cut, fun)
dem_cut2 <- calc(dem_cut, fun)




writeRasterMNT("asp_cut.asc",asp_cut2,1)
writeRasterMNT("aspect_cut.asc",aspect_cut2,1)
writeRasterMNT("slp_cut.asc",slp_cut2,1)
writeRasterMNT("lc_cut.asc",lc_cut2,1)
writeRasterMNT("st_cut.asc",st_cut2,1)
writeRasterMNT("wrf_sky_cut.asc",wrf_sky_cut2,1)
writeRasterMNT("svf_cut.asc",svf_cut2,1)
writeRasterMNT("dem_cut.asc",dem_cut2,1)



