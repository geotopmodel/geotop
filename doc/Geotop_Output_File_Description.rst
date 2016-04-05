GEOTop Output File Description
==============================

:author: Gianfranco Gallizia (eXact-lab srl)
:date: september 2014


There's a new way to define an output file in GEOtop via the
``[Output]`` section of ``geotop.inpts`` using extended keys.

Instead of specifing the file name for a pre-defined set of output files
it's now possible to define an output file for a specific variable with
a specific integration: instant value at user-defined intervals,
cumulate the values taken at user-defined intervals or average the
values on a user-defined time interval.

E.G.: Soil temperature for all layers every 24 hours::

    [Output]
    SoilTemperature::3D::INS = 86400

When GEOtop encounters those lines it will set up the output file engine
to save a group of map files (one per layer) containing the soil
temperature. The files' names will have the following format::

    201001041200_SoilTemperature_3D_INS_1d_L0001.asc
    |__________| |_____________| || | | || |   | | |
        date        variable     || | | || |   | | |
                         dimensions |_| || |   | | |
                                   type || |   | | |
                                      time |___| | |
                                          layer  |_|
                                          extension

- ``date``: the date and time in the format YYYYMMDDHHMM
- ``variable``: the name of the output variable
- ``dimensions``: one of the following
    - ``1Dp``: table with values taken at Points Of Interest (POIs)
    - ``1Ds``: table of spatial means
    - ``2D``: map
    - ``3D``: map of a single layer (part of a Tensor)
- ``type``: type of integration, one of the following
    - ``INS``: instant value
    - ``CUM``: cumulate the values taken every timestep interval and print it when you reach the right ``time``
    - ``AVG``: time average of the values
- ``time``: time of integration
- ``layer``: layer index (optional)
- ``extension``: file extension (always ``.asc``)

Setting the output directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To set the output directory for tables and maps you have to use the 
``PathPrefix`` keyword which accepts a string containing the path where 
you want to save the files.

E.G.::

    [Output]
    #Tables
    PathPrefix = "output_tables"
    SoilTemperature::1Dp::INS = 3600
    SoilTemperature::1Ds::INS = 86400
    #Maps
    PathPrefix = "output_maps"
    SoilTemperature::3D::INS = 864000

This will check if *output_tables* directory exists, if so it will be used 
for the next two files that print the POI values of temperature of the soil 
every hour and the spatial mean table every 24 hours. If *output_tables* 
directory doesn't exists GEOTop will try to create it.

The second occurence of ``PathPrefix`` sets *output_maps* as output directory 
for the next file definition that prints the soil temperature maps for each 
layer every 10 days.

Setting the layer index
^^^^^^^^^^^^^^^^^^^^^^^

To set the layer index you have to use the ``LayerIndex`` keyword.
If you set the dimension to 3D the layer index will be ignored and **all** 
the layers will be printed.

E.G.::

    [Output]
    LayerIndex = 1
    SoilTemperature::2D::INS = 86400
    SoilTemperature::3D::INS = 864000

This will print an instant value map for layer 1 every 24 hours and 
a range of maps for all layers every ten days.

The ``Layerindex`` keyword sets the value of the index only for the lines that 
**follows** it.

E.G.::

    [Output]
    SoilTemperature::3D::INS = 864000 #Valid: this will print all the layers
    SoilTemperature::2D::INS = 86400  #WARNING: this will print a map filled with NOVALUEs
    LayerIndex = 1
    SoilTemperature::2D::INS = 86400  #Valid: this will print layer 1
    LayerIndex = 2
    SoilTemperature::2D::INS = 86400  #Valid: this will print layer 2
    SoilTemperature::2D::CUM = 86400  #Valid: this will print the cumulate values of layer 2

