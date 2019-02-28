
GEOtop
######

|Build Status| |License (GPL version 3)|

:date:  last revision Fabruary 2019



**GEOtop** is a distributed model of the mass and energy balance of the
hydrological cycle, which is applicable to simulations in continuum in
small catchments. **GEOtop** deals with the effects of topography on the
interaction between energy balance and hydrological cycle with peculiar
solutions.

**GEOtop** is distributed under the GNU General Public License version 3.
A copy of the license text can be found in the COPYING file.

You can find more informations about GEOtop on the following website

                www.geotop.org 

where the model is briefly described and links to papers and other useful
websites have been collected.

Installation
************

If you want to build **GEOtop** (master branch v.2.1) from sources in your own machine:

    see here: https://github.com/geotopmodel/geotop/blob/master/doc/Install.rst 

If you prefer to install **GEOtop** via Docker to avoid manual installation of
packages:

    see here: https://hub.docker.com/r/omslab/geotop
    
    and here: https://github.com/geotopmodel/docker


Documentation
*************
    
An old version of the manual (currently under revision) can be found here:    

    http://geotopmodel.github.io/geotop/materials/geotop_manuale.pdf (updated July 2011)

in the `doc directory <https://github.com/geotopmodel/geotop/tree/master/doc>`_ there is further documentation. 
    
Documentation on former versions of the code can be found here:

    http://eprints.biblio.unitn.it/551/
    
    http://www.ing.unitn.it/dica/tools/download/Quaderni/tutorial_input_geotop.pdf
    
Useful material on **GEOtop** and his hystorical development can be found also on the R.Rigon blog:

   http://abouthydrology.blogspot.com/
   
GEOtop development branches
***************************

Currently (Agust 2018) there are several development branches in this repostory. Most used branches are the followings:

The main `**master** <https://github.com/geotopmodel/geotop>`_ branch contains the 2.1 version. It is written in c++ and it has the possibility to use the the `**MeteoIO** library <https://models.slf.ch/p/meteoio/>`_ to spatialize input meteorological variables.

This version is successfully used for operational snow mapping in the `**MySnowMaps** <http://www.mysnowmaps.com/en/>`_ app. 
However, this branch is not fully stable when the model is used with full 3D water and energy budget settings.

The `**se27xx** <https://github.com/geotopmodel/geotop/tree/se27xx>`_ branch contains the code 2.0 version which has been used for the publication  `Endrizzi et al. (2014) <https://doi.org/10.5194/gmd-7-2831-2014>`_, with some minor bug fixing. It is the most stable GEOtop version and the current benchmark for the development versions.
To install this version see https://github.com/geotopmodel/geotop/blob/se27xx/README.rst

The new `**v3.0** GEOtop development branch v3.0 (beta), written in C++, can be found in the git repo https://github.com/geotopmodel/geotop/tree/v3.0 at . You can find the compiling, running and testing instructions at https://github.com/geotopmodel/geotop/blob/v3.0/readme.md

The 3.0 version starts from version se27xx, already validated and published in the Endrizzi et al. 2014 paper.
It performs exactly as the se27xx, but it has some improvements in terms of:

- usage of object-oriented approach

- development of new data structures

- ease of compiling and running

- modularity and flexibility

- increase in testing coverage.

However, it still lacks of the integration with the MeteoIO library and other features implemented in the current 2.1 version (master branch https://github.com/geotopmodel/geotop). In the next months we plan to move toward stable 3.0 version, together with a publication. We as developer would like to have some feedbacks from you!
Bugs and suggestions can be addressed in the google groups or using pull requests.
Further test cases to validate the model are also welcome.




Report bugs/suggestion/issues
*****************************

Please use the github issues facility.

We have the following mailing lists:

   **GEOtopDev** for developers and advanced users: https://groups.google.com/forum/#!forum/geotopdev
   
   **GEOtopUsers** for regular users: https://groups.google.com/forum/#!forum/geotopusers
   

External utilities and pre-post processing scripts
**************************************************

During the years, several scripts and external softwares have been developed for preprocess **GEOtop** inputs, postprocess and visualize results. Some utilites can be found here:

There are **R scripts** (https://github.com/ecor/geotopbricks) for I/O and **GEOtop** results visualization. They work for versions 2.0 and 2.1. Mainly developed by Emanuele Cordano. There is also a stable version published on CRAN as **R package** (https://cran.r-project.org/web/packages/geotopbricks/).

The **R package Topo Sub** (https://github.com/EURAC-Ecohydro/TopoSUB) allows to produce spatially-distributed GEOtop output maps from a limited number of 1D single column simulations using a clustering approach (neglecting 3D water interactions). It has been developed by Joel Fiddes (Fiddrs and Gruber, 2012) and structured as R package by Johannes Brenner.

There are **Matlab scripts** (https://github.com/EURAC-Ecohydro/GEOmatlab) for I/O and **GEOtop** results visualization. They work for version 1.25, most of them also for version 2.0. Mainly developed by Giacomo Bertoldi and collaborators.

**GEOtop** can be embedded in the **GEOframe modelling system** (https://github.com/GEOframeOMSProjects). Mainly developed by Riccardo Rigon, Giuseppe Formetta and collaborators. For more info see: `Formetta et al. (2016a) <https://doi.org/10.3390/w8010012>`_

External models extensions 
**************************************************

Thsere are several **GEOtop** model extensions, to deal with additional physical processes. 

GEOtop model optimization tools
------------------------------

The plugin **R package geotopOtim2** (https://github.com/EURAC-Ecohydro/geotopOptim2) allows the automatic calibration and sensitivity analysis of the **GEOtop** 2.x hydrological model, based on the "Particle Swarm Optimisation" approach and the LHOAT "Latin-Hypercube One-factor-At-a-Time" approach. It has been mainly developed by Emanuele Cordano, Samuel Senoner, Giacomo Bertoldi. A paper is in preparation.

There is in development also an interface for `**PEST** <http://www.pesthomepage.org/>`_  software package for parameter estimation and uncertainty analysis. An example of the **GEOtop-PEST** interface for inverse modelling in the Rott catchment can be found at: https://doi.pangaea.de/10.1594/PANGAEA.892921. A paper by Mohsen Soltani, Patrick Laux1, Matthias Mauder, and Harald Kunstmann1 is currently in revision.

GEOtop model for shallow  landslides triggering prediction.
------------------------------
**GEOtop-SF** has been one of the first fully distributed hydrolgical models applied for hallow  landslides triggering prediction. A fundamental paper is `Simoni et al. (2008) <https://doi.org/10.1002/hyp.6886>`_, which is referred to the old 0.875 version of the model.

A more recent implementation of GEOtop for shallow landslides prectition can be found in `Formetta et al. (2016b) <https://doi.org/10.1002/2015WR017626>`_, where GEOtop is embedded in the **GEOframe modelling system**.

GEOtop model for soil erosion prediction.
------------------------------

**GEOtop_ERO** is  an extension of **GEOtop**  for modelling sediment dynamics simulating the spatio-temporal dynamics of soil erosion , deposition. Documentation can be found in `Zi et al. (2016) <https://doi.org/10.1016/j.envsoft.2016.06.004>`_

The code of the **GEOtop_ERO** model extension can be dowloaded from the repository: https://sourceforge.net/p/geotopero/code/ref/master/

GEOtop model for vegetation dynamic simulation.
------------------------------

**GEOtop_DV** is  a Matlab extension of **GEOtop**  for modelling grassland vegetation dynamics for 1D simulations. Documentation can be found in `* Della Chiesa et al. (2014) <https://doi.org/10.1002/eco.1471>`_


References
**************************************************

When using the model, please cite and refer to the following papers describing the **GEOtop** model:

* Endrizzi, S., Gruber, S., Dall’Amico, M., Rigon, R., 2014. GEOtop 2.0: simulating the combined energy and water balance at and below the land surface accounting for soil freezing, snow cover and terrain effects. Geosci. Model Dev. 7, 2831–2857. https://doi.org/10.5194/gmd-7-2831-2014

* Rigon, R., Bertoldi, G., Over, T.M., 2006. GEOtop: A Distributed Hydrological Model with Coupled Water and Energy Budgets.  J. Hydrometeorol. 7, 371–388. https://doi.org/10.1175/JHM497.1

Here is the full list of peer-reviewed publications using the GEOtop model (updated January 2019):

* Formetta, G., Capparelli, G., 2019. Quantifying the three-dimensional effects of anisotropic soil horizons on hillslope hydrology and stability. J. Hydrol. 570, 329–342. https://doi.org/10.1016/j.jhydrol.2018.12.064

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


.. |Build Status| image:: https://travis-ci.org/geotopmodel/geotop.svg?branch=master
    :target: https://travis-ci.org/geotopmodel/geotop
.. |License (GPL version 3)| image:: https://img.shields.io/badge/license-GNU%20GPL%20version%203-blue.svg
   :target: http://opensource.org/licenses/GPL-3.0



