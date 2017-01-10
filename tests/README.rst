##############
GEOtop Testing
##############


:date:  last revision january 2017 

This directory contains several tests (1D and 3D)  to check if geotop runs smoothly.
Simulation are very short and just check if code runs without any problems.

1D: contains Point simulation ( keyword PointSim=1)
3D: contains distributed simulation (keyword PointSim=0 or not defined)

Tests are automatically executed on `Travis <https://travis-ci.org/geotopmodel/geotop>`_ every push to github performed. A mail will be sent to the committer after the test is run with all the details regarding the test. You can check the log of the execution inside of Travis. 

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

Each directory does in general contains all the input files needed to run the examples.
It also contains the log files for the three referenced version.


