##############
GEOtop Testing
##############


:date:  last revision March 2017 

This directory contains several tests (1D and 3D)  to check if geotop runs smoothly.
Simulation are generally very short and just check if code runs without any problems.
T

1D: contains Point simulation ( keyword PointSim=1)
3D: contains distributed simulation (keyword PointSim=0 or not defined)

Tests are automatically executed on `Travis <https://travis-ci.org/geotopmodel/geotop>`_ every push to github performed. A mail will be sent to the committer after the test is run with all the details regarding the test. You can check the log of the execution inside of Travis. 
For each directory two kinds of test are performed: 
 a geotoop simulation is run 
 a python script is execute to compare all the output files contained in output-tabs and output-maps directories.
 

Since some test is generating an output larger than 4MB of text the complete log is disabled.

How to Execute Tests
======================

In order to perform  *output-tabs* and *output-maps* check for every test you should install **nose** and **pandas** libraries. If **pandas** is missing the files are compared exactly without considering approximation deltas.

For executing tests you have to use *ctest* command:

.. code-block:: bash

        $ ctest .

Tests are named like:

.. code-block:: text
        
        1D.<TEST_NAME>
        3D.<TEST_NAME>
        1D.<TEST_NAME>.test_runner #executes only output files check
        3D.<TEST_NAME>.test_runner #executes only output files check

You can use Regex for filtering tests and execute only what you need:

.. code-block:: bash
        
        $ ctest -R <REGEX> .

How to Add a test
=================

*CMake* is configured to check if geotop.inpts file is present in each subdirectory of *tests/1D* and *tests/3D*.
Then the test is added in ctest.
Files inside *output-tabs* and *output-maps* are compared with the outputs within the following directories, if present:

.. code-block:: text

        output-tabs-METEOIO-on
        output-maps-METEOIO-on
        output-tabs-METEOIO-off
        output-maps-METEOIO-off
        output-tabs-SE27XX
        output-maps-SE27XX

Each directory should contains all the input files needed to run the examples, a short description of the test (description.rst) and the output directory to store output. By convention these are *output-tabs* and *output-maps* .  
Once the tests is loaded into github a void placeholder files should be added in both directories. 
It also contains the log files for the three referenced version named.
A typica 1D te5st case looks like the following:

.. code-block:: text
     
      ├── DESCRIPTION.txt
      ├── geotop2-1.log-METEOIO-OFF
      ├── geotop2-1.log-METEOIO-ON
      ├── geotop.inpts
      ├── geotop.log-SE27XX
      ├── io_it.ini
      ├── listpoints.txt
      ├── meteo
      │   ├── meteo0001.txt
      ├── output-tabs
      │   ├── placeholder
      ├── output-tabs-METEOIO-OFF
      │   ├── point0001.txt
      │   ├── soilpsi0001.txt
      │   ├── soiltemp0001.txt
      │   └── thetaliq0001.txt
      ├── output-tabs-METEOIO-ON
      │   ├── point0001.txt
      │   ├── soilpsi0001.txt
      │   ├── soiltemp0001.txt
      │   └── thetaliq0001.txt
      ├── output-tabs-SE27XX
      │   ├── point0001.txt
      │   ├── soilpsi0001.txt
      │   ├── soiltemp0001.txt
      │   └── thetaliq0001.txt
      ├── rec
      ├── soil
      │ └── soil0001.txt


