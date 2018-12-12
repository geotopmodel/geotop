This directory contains  a R-shiny app to see 1D GEOtop results versus time. 
To run this app , you need to install R from CRAN website: https://www.r-project.org/ .
and install the following packages from R console previously.
Then launch R form a Shell Console: 
```
R
```
Then, once entered R environment, type the following lines to install the respective packages:

```
install.packages("devtools") ## It allows to install package directly from Github
install.packages("rgdal") 
install.packages("shiny")
install.packages("dygraphs")
install.packages("RColorBrewer")
devtools::install_github("ecor/geotopbricks") ## It is suggested to install "geotopbricks" from Github (the latest version:  1.4.9 or higher) respected to the one from CRAN!!! 


```
More details about the R packages are available on the respective CRAN page: https://cran.r-project.org/package=PACKAGE_NAME (e. g. https://cran.r-project.org/package=shiny) 
Once installed all the necessary packages, go to the app directory:
```
cd ~/[...]/tools/Test1Dshinyapp

```

Then , launch it:

```
R -e 'shiny::runApp("./",port=4566)'
```
where 4566 is the port number.
Finally, it will appear a localhost http address.
Please click on it  and view it with your own browser. 


This tool has been developed by Emanuele Cordano in 2016.
The tool has been updated by Christian Brida in 2018.

For testing the 3.0 GEOtop version versus the 2.0 version the folder should be renamed:

\Matsch_B2_Ref_007\output-tabs-SE27XX   to   \Matsch_B2_Ref_007\output-tabs-v_2.0

\Matsch_B2_Ref_007\output-tabs-METEOIO-OFF   to   \Matsch_B2_Ref_007\output-tabs-v_3.0


