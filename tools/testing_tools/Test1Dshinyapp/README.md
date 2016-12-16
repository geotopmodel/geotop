This directory contains several an R-shiny app to see 1D GEOtop results versus time. 
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
devtools::install_github("ecor/geotopbricks") ## It is suggested to install "geotopbricks" from Github (the latest version) respected to the one from CRAN!!! 


```
More details about the R packages are available on the respective CRAN page: https://cran.r-project.org/package=PACKAGE_NAME (e. g. https://cran.r-project.org/package=shiny) 
Once installed all the necessary packages, go to the app directory:
```
cd ~/[...]/tools/Test1Dshinyapp

```

Then , launch it:

```
R -e 'shiny::runApp("./")'
```

Finally, it will appear a localhost http address.
Please click on it  and view it with your own browser. 
