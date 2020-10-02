#  GEOtop version geotop20_veg

This is a branch of GEOtop 2.0 with some bug fixes and added input and output features to improve vegetation modelling.
The purpose of this branch is to improve evapotranspiration parametrization for GEOtop,
with particular reference to mountain grassland  (LTER and Monalisa test sites managed by EURAC)
and apple orchards (Caldaro test site managed by LUB).

Prepared by G. Bertoldi and S. Senoner

In this file are described the major changes with respect GEOtop 2.0 (branch se27xx)

## How to compile

type from this path:

mkdir bin # This is the directory for the executable file
make -f geotop.make

This will put the executable in the subdirectory "bin"

## Major updates

### added averaged SWE to basin output 

https://github.com/geotopmodel/geotop/commit/326943d165b59a8d52eb47a25f622728acf1d727

### Introduced new keywords for vegetation transpiration parameters controlling vegetation stress

### keywords as flag to  control Jarvis type stomatal representation 

VegRswStress  //(default =1, solar radiation stress [Best, (1998); Dolman et al., 1991])
VegVPDStress //(default 1= air water vapour pressure deficit stress [Best, (1998); Dickinson et al., 1991])
VegTempStress // (default =1 air temperature stress [Best, (1998); Dickinson et al., 1991])
VegWaterStress // (default =1 soil water content stress [Wigmosta et al., (1994); Feddes et al.(1978)])
RsLAI // (default =0 ) if(par->RsLAI == 1){ Rsmin=land[jrs]/LSAI; } else { Rsmin=land[jrs];		}

### keywords for parameters to control vegetation stomata VPD stress following Dickinson et al., 1991

VpdvegMax
// Vegetation Vapor Pressure Deficit  stomata stress factor (default 40[hPa]) fe=1.0-(ev-ea)/VpdvegMax

### keywords for parameters to control vegetation stomata temperature stress 
following Dickinson et al., 1991 fTemp=(Tv-TvegMin_)*( TvegMax-Tv)/TvegRes;

TvegMin  // Minumum working leaves temperature for stomata default 0 [C]
TvegMax // Maximum working leaves temperature for stomata default 50 [C]
TvegRes //Stomata temperature stress factor default 625  [C^2]

## Bibliography

Feddes, R. A., Kowalik, P. J., & Zaradny, H. (1978). Simulation of field water use and crop yield. In Simulation Monographs (p. 188pp.). Wageningen: PUDOC.

Jarvis, P. G., & Morrison, J. I. L. (1981). The control of transpiration and photosynthesis by the stomata. In P. G. Jarvis & T. A. Mansfield (Eds.), Stomatal Physiology (pp. 247-279). UK: Cambridge Univ. Press.

Best, M. J. (1998). A model to Predict Surface Temperatures. Bound. Layer Meteorol., 88(2), 279–306.

Dickinson, R.E., A. Henderson-Sellers, C. Rosenzweig, and P.J. Sellers, 1991: Evapotranspiration models with canopy resistance for use in climate models, a review. Agric. Forest Meteorol., 54, 373-388, doi:10.1016/0168-1923(91)90014-H.

Dolman, A. J., Gash, J. H. C., Roberts, J., & Shuttleworth, W. J. (1991). Stomatal and surface conductance of tropical rain-forest. Agricultural and Forest Meteorology, 54, 303−318

Wigmosta, M. S., Vail, L., & Lettenmaier, D. (1994). A Distributed Hydrology-Vegetation Model for complex terrain. Water Resour. Res., 30(6), 1665–1679.
