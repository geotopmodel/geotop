# GEOtop 3.0
######

|Build Status| |License (GPL version 3)|

:date:  last revision April 2022


***********************************************************
## Introduction

**GEOtop** is a distributed model of the mass and energy balance of the
hydrological cycle, which is applicable to simulations in continuum in
small catchments. **GEOtop** deals with the effects of topography on the
interaction between energy balance and hydrological cycle with peculiar
solutions.

**GEOtop** is distributed under the GNU General Public License version 3.
A copy of the license text can be found in the COPYING file.


**GEOtop 3.0** version (www.geotop.org) starts from 2.0 version (branch ```se27xx```) already validated and published
in the Endrizzi et al. 2014 paper.
It performs exactly as the previous version but it has some improvements in terms of:
- usage of object-oriented approach
- development of new data structures
- ease of compiling and running
- modularity and flexibility
- increase in testing coverage.

However, the 3.0 version  still lacks of the integration with the MeteoIO library and other features implemented in the current [2.1.1 version](https://github.com/geotopmodel/geotop/tree/2.1.1) (former branch ```master```).
In the next months we plan to move toward a stable 3.0 version, together with a publication.

A more detailed description of this new version can be found in the MHPC thesis (https://www.mhpc.it/)
of Elisa Bortoli at the following link: https://iris.sissa.it/handle/20.500.11767/86154#.XEXzOsZ7l8w

You can compile and run GEOtop 3.0 following the listed instructions,
tested on a Linux system.

***********************************************************
## Getting the source code

- Clone the git repository:
```
git clone https://github.com/geotopmodel/geotop.git
```

- Move to your local geotop repository:
```
cd geotop
```

- Go into the branch v3.0 and make sure you are in:
```
git checkout v3.0
git branch
```

- to update your local repository to the newest commit, execute 
```
git pull origin v3.0
```

***********************************************************
## Compiling
Now you can compile using a build system tool.
Build tools are programs that automate the creation of executable applications
from source code.
Building incorporates compiling, linking and packaging the code into
a usable or executable form.
You can choose between [CMake](https://cmake.org/) (Please note that you need cmake version 3 at least)
and [Meson](http://mesonbuild.com/).

### Option 1: using CMake

Note that you need to have installed in your system CMake version > 3.0 as it's written in the CMakeLists.txt (https://github.com/geotopmodel/geotop/blob/v3.0/CMakeLists.txt).

- Create the build directory and go inside it:
```
mkdir cmake-build
cd cmake-build
```
- Check the default values for the options, opening the file CMakeList.txt
in the upper directory or writing:
```
ccmake ..
```
- Press [c] and [e] to configure and edit the options

- Press [t] to toggle the advanced mode; several options will appear.
You can modify the values of the flags going to the correspondent line,
pressing "Enter" key and then editing; to save what you have just written
press again "Enter".

- For example you can choose the build type, writing RELEASE (default option) or
DEBUG after ```CMAKE_BUILD_TYPE```,
and modify the other flags as you prefer, knowing that a flag like:
    - *_RELEASE: will be applied only when compiling in RELEASE mode
    - *_DEBUG: will be applied only when compiling in DEBUG mode.

- Press again [c] and [e] to configure; then press [g] to generate and exit.
Now the current directory will have the following files and folders:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/make-build[v3.0*] $ ls -l
totale 84
-rw-rw-r-- 1 elisa elisa 11518 giu 15 14:22 CMakeCache.txt
drwxrwxr-x 5 elisa elisa  4096 giu 15 14:23 CMakeFiles
-rw-rw-r-- 1 elisa elisa  1574 giu 15 14:23 cmake_install.cmake
-rw-rw-r-- 1 elisa elisa   307 giu 15 14:23 CTestTestfile.cmake
-rw-rw-r-- 1 elisa elisa 51860 giu 15 14:23 Makefile
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:22 src
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:23 tests
```
- Compile (-j4 allows the usage of 4 processes):
```
make -j4
```

### Option 2: using Meson
- Create the build directory and go inside it:
```
mkdir meson-build
cd meson-build
```

- Create the build file:
```
meson
```

- Now the current directory will have the following files and folders:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build[v3.0*] $ ls -l
totale 156
-rw-rw-r-- 1 elisa elisa 64560 giu 15 14:39 build.ninja
-rw-rw-r-- 1 elisa elisa 67030 giu 15 14:39 compile_commands.json
drwxrwxr-x 7 elisa elisa  4096 giu 15 14:39 meson
drwxrwxr-x 2 elisa elisa  4096 giu 15 14:39 meson-logs
drwxrwxr-x 2 elisa elisa  4096 giu 15 14:39 meson-private
drwxrwxr-x 4 elisa elisa  4096 giu 15 14:39 src
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:39 subprojects
drwxrwxr-x 3 elisa elisa  4096 giu 15 14:39 tests
```

- Check the default values for the options, opening the file meson.build
in the upper directory or typing:
```
meson configure
```

- If you want to modify some of them, add -Doption=value: for example
    - to set the build type to debug writes:
    ```
    meson configure -Dbuildtype=debug
    ```
    - to add compiler and linker options (i.e. add ```-pg```) write:
    ```
    meson configure -Dcpp_args=-pg -Dcpp_link_args=-pg
    ```
    - define multiple compiler options:
    ```
    meson configure -Dcpp_args=" -OPTION_1 -OPTION_2"
    ```

- To check if the desired flags were activated, you can look at their current values
(```true``` or ```false```) again typing inside the build folder:
```
meson configure
```

- Compile:
```
ninja
```

***********************************************************
## Running GEOtop

If you have compiled the code with CMake you will have your executable called ```geotop``` in the folder ```cmake-build```
with meson in the folder ```meson-build```.

To run the code just type

```
./geotop simulation_path
```
where simulation_path if the path to a folder where are all the input data of a simulation for example:

```
./geotop ../tests/1D/Matsch_B2_Ref_007/
```

You can then create a link to the executable or install itclsls

***********************************************************
## Running the test cases
Now you can run the proposed test cases.

Note that you need to have **numdiff** installed on your mashine 
```
apt-get install -y numdiff
```

### Option 1: using CMake
- Know which tests are available:
```
 ctest -N
```

- Run a single test (i.e. Mazia):
```
ctest -R Mazia
```

- Run a group of tests (i.e. all 1D tests, using 4 processes):
```
ctest -R 1D -j4
```

- Run all tests
```
ctest
```



### Option 2: using Meson
- Know which tests are available:
```
 meson test --list
 ```

- Run a single test (i.e. Mazia):
```
meson test --suite geotop:Mazia
```

- Run a group of tests (i.e. all 1D tests):
```
meson test --suite geotop:1D
```

- Run all tests
```
ninja test
```

***********************************************************
## Getting help

### Output error message
If for some reasons at a certain point after typing ```ninja``` you get a message like:
```
Something went terribly wrong. Please file a bug.
FAILED: build.ninja
```
remove the build folder and create it again.

### Looking into the documentation
An interactive documentation can be built with Doxygen
by typing in the root directory:
```
doxygen Doxyfile
```
A new folder **doxygen_generated_doc** will be created, containing two subfolders:
*html* and *latex*, with graphs of the single functions.

If you want to navigate files and functions, go inside
*html* folder and type:
```
firefox index.html
```

### Reporting an issue
To report a problem you can open an issue on GitHub (https://github.com/geotopmodel/geotop/issues) listing all the following infos:
- short description of what happened
- operating system (OS)
- architecture
- compiler
- configuration (type from the build folder ```ccmake ..```
  or ```meson configure```, depending on the build system tool you are using,
  and put the info in a file).
  
### GEOtop users community

If you want to get in contact with the users and developers community or discuss about your GEOtop application we have the following **mailing lists**:

   **GEOtopDev** for developers and advanced users: https://groups.google.com/forum/#!forum/geotopdev
   
   **GEOtopUsers** for regular users: https://groups.google.com/forum/#!forum/geotopusers
  
***********************************************************
## DOCUMENTATION
  
    
A  **manual** of the model (for the version 1.2, mainly valid also for the current version) can be found here:    

http://geotopmodel.github.io/geotop/materials/geotop_manuale.pdf (updated July 2011)

in the [doc directory](https://github.com/geotopmodel/geotop/tree/master/doc) there is further documentation. 
    
Documentation on former versions of the code can be found here:

http://eprints.biblio.unitn.it/551/
    
http://www.ing.unitn.it/dica/tools/download/Quaderni/tutorial_input_geotop.pdf
    
Useful material on **GEOtop** and his hystorical development can be found also on the R.Rigon blog:

http://abouthydrology.blogspot.com/

*********************************************************
## GEOTOP MODEL I/0 SCRIPTING TOOLS and PLUGINS

The GEOtop model used as I/0 ascii text files. To better analyze and exploit GEOtop outputs several tools have been prepared.

### **R inteface** with the geotopbricks R package 

The [**geotopbricks (R Package)**](https://CRAN.R-project.org/package=geotopbricks) develope by E. Cordano allow a full set of functions to integrate GEOtop outptus as data strctures in the R environment.
The development version is in the following git repository https://github.com/ecor/geotopbricks

### **Python wrapper** with the GEOtoPy Python package

[**GEOtoPy**](https://github.com/stefanocampanella/GEOtoPy) is a small Python package developed by S. Campanella. It works as a GEOtop Python wrapper, exporting a single base class GEOtoPy.GEOtop.

### **Matlab scripts** with the GEOmatlab collection

[**GEOmatlab**](https://github.com/EURAC-Ecohydro/GEOmatlab) is a collection of Matlab scripts to import, analyze and plot GEOtop model output files developed by G. Bertoldi. 
   
*********************************************************
## EXTERNAL MODEL EXTENSIONS 

Thsere are several **GEOtop** model extensions, to deal with additional physical processes. 

### **High-Performance Optimization** for the Calibration of the GEOtop Model

The repository [**Stefanocampanella/MHPC-project**](https://github.com/stefanocampanella/MHPC-project) contains notebooks, code and documentation for a high-performance derivative-free optimization to exploit HPC  for the calibration of parameters of the **GEOtop** model. It has been developed by Stefano Campanella in the course of his MHPC Thesis [**Calibration of the GEOtop model using evolutionary algorithms on supercomputers**](https://stefanocampanella.github.io/MHPC-project/home.html)

### GEOtop model particle swarm optimization with **R**

The plugin [**geotopOtim2 (R Package)**](https://github.com/EURAC-Ecohydro/geotopOptim2), based on [**geotopbricks (R Package)**](https://CRAN.R-project.org/package=geotopbricks) allows the automatic calibration and sensitivity analysis of the **GEOtop** 2.x hydrological model, based on the "Particle Swarm Optimisation" approach and the LHOAT "Latin-Hypercube One-factor-At-a-Time" approach. It has been mainly developed by Emanuele Cordano, Samuel Senoner, Giacomo Bertoldi. 

### GEOtop model optimization with **PEST**  

It has been developed an  interface for  [**PEST**](http://www.pesthomepage.org/)  software package for parameter estimation and uncertainty analysis. An example of the **GEOtop-PEST** interface for inverse modelling in the Rott catchment can be found at: https://doi.pangaea.de/10.1594/PANGAEA.892921. Full details can be found in the paper  [Soltani et al. (2019)](https://doi.org/10.1016/j.jhydrol.2019.02.033)

In general, PEST requires the following input files for automatic parameter estimation and inverse modelling: (i) Template files, to identify the model parameters; (ii) Instruction files, to identify the model outputs; and (iii) Control file, which supplies PEST with the names of all template and instruction files, the names of model input and output files, initial parameter values, measurement values and weights, etc. (Doherty, 2010). 

The  PEST software (Doherty, 2002) together with over 100-utility-programs such as SENSAN
and GENLINPRED used herein are freely available at  http://www.pesthomepage.org/Downloads.php.  For detailed and comprehensive information for combining a model of interest with PEST, it is referred to Sect. “3. The Model-PEST Interface” of the PEST manual, as described in Doherty (2002).

### GEOtop model for shallow  landslides triggering prediction.

**GEOtop-SF** has been one of the first fully distributed hydrolgical models applied for hallow  landslides triggering prediction. A fundamental paper is [Simoni et al. (2008)](https://doi.org/10.1002/hyp.6886), which is referred to the old 0.875 version of the model.

A more recent implementation of GEOtop for shallow landslides prectition can be found in [Formetta et al. (2016b)](https://doi.org/10.1002/2015WR017626), where GEOtop is embedded in the **GEOframe modelling system**.

### GEOtop model for soil erosion prediction.

**GEOtop_SED** is  an extension of **GEOtop**  for modelling sediment dynamics simulating the spatio-temporal dynamics of soil erosion , deposition. Documentation can be found in [Zi et al. (2016)](https://doi.org/10.1016/j.envsoft.2016.06.004)

The code of the **GEOtop_sed** model extension can be dowloaded from the repository: 
https://github.com/TanZiTT/GEOtopSed

### GEOtop model for vegetation dynamic simulation.

**GEOtop_DV** is  a Matlab extension of **GEOtop**  for modelling grassland vegetation dynamics for 1D simulations. Documentation can be found in [Della Chiesa et al. (2014)](https://doi.org/10.1002/eco.1471)

**************************************************
## Operational **GEOtop** model applications 

The GEOtop model has been also used for operational application:

### Snow depth mapping

The **GEOtop** model (v 2.1) is the scientific basis of the [**MySnowMaps**](https://www.mysnowmaps.com/en/) service, which presents real time snow depth maps and prediction for the Alps, implemented by M. Dall´Amico the [**MobyGis**](http://www.mobygis.com/)  company.

The **GEOtop** model (v 2.1) has been also used by P. Pogliotti for the ARPA Valle d´Aosta (Italy) to monitor in real-time meltwater avaliabilty for hydropower.


### Water budget mapping

A preliminary application of the **GEOtop** model (v 3.0) for mapping the water budget of the Venosta (Italy) catchment in near real time on a weekly basis has implemented in the following web-gis: https://maps.civis.bz.it/ in the framework of the European Regional Development Fund (ERDF) project DPS4ESLAB.


***********************************************************
## ACKNOWLEDGEMENTS

### Mandatory references

When using the model, **please cite the following fundamental papers** describing the **GEOtop** model:

* Endrizzi, S., Gruber, S., Dall’Amico, M., Rigon, R., 2014. GEOtop 2.0: simulating the combined energy and water balance at and below the land surface accounting for soil freezing, snow cover and terrain effects. Geosci. Model Dev. 7, 2831–2857. https://doi.org/10.5194/gmd-7-2831-2014

* Rigon, R., Bertoldi, G., Over, T.M., 2006. GEOtop: A Distributed Hydrological Model with Coupled Water and Energy Budgets.  J. Hydrometeorol. 7, 371–388. https://doi.org/10.1175/JHM497.1

### GEOtop developers

The GEOtop model has been developed since year 2000 by a number of people, starting from the research group of Prof. R. Rigon of the [University of Trento, Italy](https://www.dicam.unitn.it/en), and then by different reseach groups worldwide, in particular the [University of Zurich, CH](https://www.uzh.ch/en.html), [Eurac research, Italy](www.eurac.edu), the Mountaneering, MobyGIS , [WaterJade](https://www.waterjade.com/en/mobygis-en/), [Rendena100](https://www.rendena100.eu) )companies.   A non exaustive list of contributors include: Marco Pegoretti, Giacomo Bertoldi, Fabrizio Zanotti, Silvia Simoni, Stefano Endrizzi, Matteo dell´Amico, Emanuele Cordano , Stefan Gruber, Andrea Cozzini, Alberto Sartori, Samuel Senoner, Elisa Bortoli.

### Financial support

Part of the GEOtop v3.0 refactoring work was developed by E. Bortoli in the framework of her MHCP thesis.
This last work was partly supported by:
- OGS, CINECA and EURAC Research under HPC-TRES program award number 2017-20.
- The European Regional Development Fund, Operational Programme Investment for growth and jobs ERDF 2014-2020 under Project number ERDF1094, Data Platform and Sensing Technology for Environmental Sensing LAB – DPS4ESLAB.”

**************************************
## REFERENCES

Here is the full list of peer-reviewed publications using the GEOtop model (updated October 2021):

* Wani, J. M., Thayyen, R. J., Ojha, C. S. P., and Gruber, S.: The surface energy balance in a cold and arid permafrost environment, Ladakh,  Himalayas, India, The Cryosphere, 15, 2273--2293, https://doi.org/10.5194/tc-15-2273-2021, 2021.

* Bright Ross, J.G., Peters, W., Ossi, F., Moorcroft, P.L.,  Cordano, E.,  Eccel, E.,  Bianchini, F.,  Ramanzin, M., and  Cagnacci, F. . Climate change and anthropogenic food manipulation interact in shifting the distribution of a large herbivore at its altitudinal range limit. Sci Rep 11, 7600 (2021). https://doi.org/10.1038/s41598-021-86720-2

*  Terzago, S., Andreoli, V., Arduini, G., Balsamo, G., Campo, L., Cassardo, C., Cremonese, E., Dolia, D., Gabellani, S., Hardenberg, J. von, Cella, U.M. di, Palazzi, E., Piazzi, G., Pogliotti, P., Provenzale, A., 2020. Sensitivity of snow models to the accuracy of meteorological forcings in mountain environments. Hydrology and Earth System Sciences 24, 4061–4090. https://doi.org/10.5194/hess-24-4061-2020

* Wani, J.M., Thayyen, R.J., Gruber, S., Ojha, C.S.P., Stumm, D., 2020. Single-year thermal regime and inferred permafrost occurrence in the upper Ganglass catchment of the cold-arid Himalaya, Ladakh, India. Sci. Total Environ. 703, 134631. https://doi.org/10.1016/j.scitotenv.2019.134631

* Zi, T., Kumar, M., Albertson, J., 2019. Intercomparing varied erosion, deposition and transport process representations for simulating sediment yield. Sci. Rep. 9, 1–13. https://doi.org/10.1038/s41598-019-48405-9

* Fiddes, J., Aalstad, K., Westermann, S., 2019. Hyper-resolution ensemble-based snow reanalysis in mountain regions using clustering. Hydrol. Earth Syst. Sci. 23, 4717–4736. https://doi.org/10.5194/hess-23-4717-2019

* Fullhart, A.T., Kelleners, T.J., Speckman, H.N., Beverly, D., Ewers, B.E., Frank, J.M., Massman, W.J., 2019. Measured and Modeled Above‐ and Below‐Canopy Turbulent Fluxes for a Snow‐Dominated Mountain Forest Using Geotop, Hydrological Processes. https://doi.org/10.1002/hyp.13487

* Soltani, M., Laux, P., Mauder, M., Kunstmann, H., 2019. Inverse distributed modelling of streamflow and turbulent fluxes: A sensitivity and uncertainty analysis coupled with automatic optimization. J. Hydrol. 571, 856–872. https://doi.org/10.1016/j.jhydrol.2019.02.033

* Formetta, G., Capparelli, G., 2019. Quantifying the three-dimensional effects of anisotropic soil horizons on hillslope hydrology and stability. J. Hydrol. 570, 329–342. https://doi.org/10.1016/j.jhydrol.2018.12.064

* Kiese, R., Fersch, B., Baessler, C., Brosy, C., Butterbach-Bahl, K., Chwala, C., Dannenmann, M., Fu, J., Gasche, R., Grote, R., Jahn, C., Klatt, J., Kunstmann, H., Mauder, M., Rödiger, T., Smiatek, G., Soltani, M., Steinbrecher, R., Völksch, I., Werhahn, J., Wolf, B., Zeeman, M., Schmid, H.P., 2018. The TERENO Pre-Alpine Observatory: Integrating Meteorological, Hydrological, and Biogeochemical Measurements and Modeling. Vadose Zo. J. 17, 0. https://doi.org/10.2136/vzj2018.03.0060

* Soltani, M., Laux, P., Mauder, M., Kunstmann, H., 2018. Spatiotemporal variability and empirical Copula-based dependence structure of modeled and observed coupled water and energy fluxes. Hydrol. Res. nh2018163. https://doi.org/10.2166/nh.2018.163

* Pullens, J.W.M., Sottocornola, M., Kiely, G., Gianelle, D., Rigon, R., 2018. Assessment of the water and energy budget in a peatland catchment of the Alps using the process based GEOtop hydrological model. J. Hydrol. 563, 195–210. https://doi.org/10.1016/j.jhydrol.2018.05.041

* Fullhart, A.T., Kelleners, T.J., Chandler, D.G., Mcnamara, J.P., Seyfried, M.S., 2018. Water Flow Modeling with Dry Bulk Density Optimization to Determine Hydraulic Properties in Mountain Soils. Soil Sci. Soc. Am. J. 82, 31–44. https://doi.org/10.2136/sssaj2017.06.0196

* Kollet, S., Sulis, M., Maxwell, R.M.R.M., Paniconi, C., Putti, M., Bertoldi, G., Coon, E.T.E.T., Cordano, E., Endrizzi, S., Kikinzon, E., Mouche, E., Mügler, C., Park, Y.-J.Y.-J., Refsgaard, J.C.J.C., Stisen, S., Sudicky, E., 2017. The integrated hydrologicmodel intercomparison project, IH-MIP2: A second set of benchmark results to diagnose integrated hydrology and feedbacks. Water Resour. Res. 53, 867–890. https://doi.org/10.1002/2014WR015716

* Engel, M., Notarnicola, C., Endrizzi, S., Bertoldi, G., 2017. Snow model sensitivity analysis to understand spatial and temporal snow dynamics in a high-elevation catchment. Hydrol. Process. 31, 4151–4168. https://doi.org/10.1002/hyp.11314

* Mauder, M., Genzel, S., Fu, J., Kiese, R., Soltani, M., Steinbrecher, R., Zeeman, M., Banerjee, T., De Roo, F., Kunstmann, H., 2017. Evaluation of energy balance closure adjustment methods by independent evapotranspiration estimates from lysimeters and hydrological simulations. Hydrol. Process. https://doi.org/10.1002/hyp.11397

* Formetta, G., Capparelli, G., David, O., Green, T.R., Rigon, R., 2016. Integration of a Three-Dimensional Process-Based Hydrological Model into the Object Modeling System. Water 8, 1–15. https://doi.org/10.3390/w8010012

* Hingerl, L., Kunstmann, H., Wagner, S., Mauder, M., Bliefernicht, J., Rigon, R., 2016. Spatio-temporal variability of water and energy fluxes - a case study for a mesoscale catchment in pre-alpine environment. Hydrol. Process. 30, 3804–3823. https://doi.org/10.1002/hyp.10893

* Zi, T., Kumar, M., Kiely, G., Lewis, C., Albertson, J., 2016. Simulating the spatio-temporal dynamics of soil erosion , deposition , and yield using a coupled sediment dynamics and 3D distributed hydrologic model. Environ. Model. Softw. 83, 310–325. https://doi.org/10.1016/j.envsoft.2016.06.004

* Formetta, G., Simoni, S., Godt, J.W., Lu, N., Rigon, R., 2016. Geomorphological control on variably saturated hillslope hydrology and slope instability. Water Resour. Res. 52, 4590–4607. https://doi.org/10.1002/2015WR017626

* Greifeneder, F., Notarnicola, C., Bertoldi, G., Brenner, J., Wagner, W., 2015. A novel approach to improve spatial detail in modeled soil moisture through the integration of remote sensing data, in: Geoscience and Remote Sensing Symposium (IGARSS), 2015 IEEE International. pp. 1988–1991. https://doi.org/10.1109/IGARSS.2015.7326187

* Fiddes, J., Endrizzi, S., Gruber, S., 2015. Large-area land surface simulations in heterogeneous terrain driven by global data sets : application to mountain permafrost. Cryosph. 9, 411–426. https://doi.org/10.5194/tc-9-411-2015

* Eccel, E., Cordano, E., Zottele, F., 2015. A project for climatologic mapping of soil water content in Trentino. Ital. J. Agrometeorol. 1, 5–20.

* Bertoldi, G., Della Chiesa, S., Notarnicola, C., Pasolli, L., Niedrist, G., Tappeiner, U., Della, S., Notarnicola, C., Pasolli, L., Niedrist, G., Tappeiner, U., 2014. Estimation of soil moisture patterns in mountain grasslands by means of SAR RADARSAT2 images and hydrological modeling. J. Hydrol. 516, 245–257. https://doi.org/10.1016/j.jhydrol.2014.02.018

* Endrizzi, S., Gruber, S., Dall’Amico, M., Rigon, R., 2014. GEOtop 2.0: simulating the combined energy and water balance at and below the land surface accounting for soil freezing, snow cover and terrain effects. Geosci. Model Dev. 7, 2831–2857. https://doi.org/10.5194/gmd-7-2831-2014

* Della Chiesa, S., Bertoldi, G., Niedrist, G., Obojes, N., Endrizzi, S., Albertson, J.D., Wohlfahrt, G., Hörtnagl, L., Tappeiner, U., 2014. Modelling changes in grassland hydrological cycling along an elevational gradient in the Alps. Ecohydrology n/a--n/a. https://doi.org/10.1002/eco.1471

* Cordano, E., Rigon, R., 2013. A mass-conservative method for the integration of the two-dimensional groundwater (Boussinesq) equation. Water Resour. Res. 49, 1058–1078. https://doi.org/10.1002/wrcr.20072

* Lewis, C., Albertson, J., Zi, T., Xu, X., Kiely, G., 2013. How does afforestation affect the hydrology of a blanket peatland? A modelling study. Hydrol. Process. 27, 3577–3588. https://doi.org/10.1002/hyp.9486

* Gubler, S., Endrizzi, S., Gruber, S., Purves, R.S., 2013. Sensitivities and uncertainties of modeled ground temperatures in mountain environments. Geosci. Model Dev. 6, 1319–1336. https://doi.org/10.5194/gmd-6-1319-2013

* Fiddes, J., Gruber, S., 2012. TopoSUB: a tool for efficient large area numerical modelling in complex topography at sub-grid scales. Geosci. Model Dev. 5, 1245–1257. https://doi.org/10.5194/gmd-5-1245-2012

* Dall’Amico, M., Endrizzi, S., Gruber, S., Rigon, R., 2011. A robust and energy-conserving model of freezing variably-saturated soil. Cryosph. 5, 469–484. https://doi.org/10.5194/tc-5-469-2011

* Bertoldi, G., Notarnicola, C., Leitinger, G., Endrizzi, S., Della Chiesa, S., Zebisch, M., Tappeiner, U., Della Chiesa, S., Tappeiner, U., 2010. Topographical and ecohydrological controls on land surface temperature in an Alpine catchment. Ecohydrology 3, 189–204. https://doi.org/10.1002/eco.129

* Endrizzi, S., Marsh, P., 2010. Observations and modeling of turbulent fluxes during melt at the shrub-tundra transition zone 1: point scale variations. Hydrol. Res. 41, 471–490.

* Gebremichael, M., Rigon, R., Bertoldi, G., Over, T.M.M., 2009. On the scaling characteristics of observed and simulated spatial soil moisture fields. Nonlin. Process. Geophys. 16, 141–150. https://doi.org/10.5194/npg-16-141-2009

* Simoni, S., Zanotti, F., Bertoldi, G., Rigon, R., 2008. Modelling the probability of occurrence of shallow landslides and channelized debris flows using GEOtop-FS. Hydrol. Process. doi: 10.10, 532–545. https://doi.org/10.1002/hyp.6886

* Bertoldi, G., Rigon, R., Over, T.M.M., 2006. Impact of Watershed Geomorphic Characteristics on the Energy and Water Budgets. J. Hydrometeorol. 7, 389–403. https://doi.org/10.1175/JHM500.1

* Rigon, R., Bertoldi, G., Over, T.M.M., 2006. GEOtop: A Distributed Hydrological Model with Coupled Water and Energy Budgets. J. Hydrometeorol. 7, 371–388. https://doi.org/10.1175/JHM497.1

* Zanotti, F., Endrizzi, S., Bertoldi, G., Rigon, R., 2004. The GEOtop snow module. Hydrol. Proc. 18, 3667–3679. DOI:10.1002/hyp.5794. https://doi.org/10.1002/hyp.5794


