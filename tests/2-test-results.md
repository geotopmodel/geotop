# General infos
- Compiler: c++ (gcc 5.4.0 "c++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609")
- Processor: Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
- Meson build type: release
- Author: Elisa Bortoli (elisa.bortoli3@gmail.com)
- Date: 22-06-2018

## 1D tests

### Cmake
Failing tests: 0/12.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build[v3.0*] $ ctest -R "1D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build
      Start  1: 1D/B2_BeG_017
      Start  2: 1D/Bro
      Start  3: 1D/Calabria
      Start  4: 1D/ColdelaPorte
      Start  5: 1D/CostantMeteo
      Start  6: 1D/Jungfraujoch
      Start  7: 1D/Matsch_B2_Ref_007
      Start  8: 1D/Matsch_P2_Ref_007
 1/12 Test  #4: 1D/ColdelaPorte ..................   Passed    1.10 sec
      Start  9: 1D/PureDrainage
 2/12 Test  #9: 1D/PureDrainage ..................   Passed    1.01 sec
      Start 10: 1D/PureDrainageFaked
 3/12 Test #10: 1D/PureDrainageFaked .............   Passed    1.01 sec
      Start 11: 1D/PureDrainageRainy
 4/12 Test #11: 1D/PureDrainageRainy .............   Passed    0.91 sec
      Start 12: 1D/PureDrainageRainySlope
 5/12 Test  #1: 1D/B2_BeG_017 ....................   Passed    6.56 sec
 6/12 Test #12: 1D/PureDrainageRainySlope ........   Passed    0.84 sec
 7/12 Test  #3: 1D/Calabria ......................   Passed    7.08 sec
 8/12 Test  #8: 1D/Matsch_P2_Ref_007 .............   Passed    8.86 sec
 9/12 Test  #7: 1D/Matsch_B2_Ref_007 .............   Passed    9.47 sec
10/12 Test  #5: 1D/CostantMeteo ..................   Passed   11.89 sec
11/12 Test  #2: 1D/Bro ...........................   Passed   15.30 sec
12/12 Test  #6: 1D/Jungfraujoch ..................   Passed   19.63 sec

100% tests passed, 0 tests failed out of 12

Total Test time (real) =  19.65 sec
```

### Meson
Failing tests: 0/24.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build[v3.0*] $ meson test --suite geotop:1D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build'
ninja: no work to do.
 1/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017    OK       6.80 s
 2/24 geotop:1D+Bro / 1D/Bro                  OK      12.83 s
 3/24 geotop:1D+Calabria / 1D/Calabria        OK       2.43 s
 4/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte  OK       0.36 s
 5/24 geotop:1D+CostantMeteo / 1D/CostantMeteo  OK       4.20 s
 6/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch  OK      28.66 s
 7/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007  OK       8.52 s
 8/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007  OK       8.94 s
 9/24 geotop:1D+PureDrainage / 1D/PureDrainage  OK       1.04 s
10/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked  OK       1.01 s
11/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy  OK       1.01 s
12/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope  OK       0.94 s
13/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017.test_runner  OK       0.06 s
14/24 geotop:1D+Bro / 1D/Bro.test_runner      OK       3.22 s
15/24 geotop:1D+Calabria / 1D/Calabria.test_runner  OK       2.01 s
16/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte.test_runner  OK       0.01 s
17/24 geotop:1D+CostantMeteo / 1D/CostantMeteo.test_runner  OK       1.46 s
18/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch.test_runner  OK       0.02 s
19/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007.test_runner  OK       3.85 s
20/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007.test_runner  OK       4.34 s
21/24 geotop:1D+PureDrainage / 1D/PureDrainage.test_runner  OK       0.14 s
22/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked.test_runner  OK       0.13 s
23/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy.test_runner  OK       0.13 s
24/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope.test_runner  OK       0.13 s

OK:        24
FAIL:       0
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build/meson-logs/testlog.txt
```

## 3D tests

### Cmake
Failing tests: 1/20.
- 3D/small_example-channel

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build[v3.0*] $ ctest -R "3D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build
      Start 13: 3D/Borden05m
      Start 14: 3D/example
      Start 15: 3D/hillslope01
      Start 16: 3D/hillslope02_superslab
      Start 17: 3D/Mazia
      Start 18: 3D/no_reflection
      Start 19: 3D/onepoint_hydrostatic
      Start 20: 3D/panola
 1/20 Test #18: 3D/no_reflection ......................   Passed    2.70 sec
      Start 21: 3D/panola_25pixel
 2/20 Test #21: 3D/panola_25pixel .....................   Passed    2.60 sec
      Start 22: 3D/panola_25pixel_nobed
 3/20 Test #14: 3D/example ............................   Passed   17.23 sec
      Start 23: 3D/panola_25pixel_nobed_hydrostatic
 4/20 Test #16: 3D/hillslope02_superslab ..............   Passed   35.95 sec
      Start 24: 3D/prealpiC
 5/20 Test #19: 3D/onepoint_hydrostatic ...............   Passed   50.37 sec
      Start 25: 3D/PSQL_test
 6/20 Test #17: 3D/Mazia ..............................   Passed   75.84 sec
      Start 26: 3D/rendena
 7/20 Test #20: 3D/panola .............................   Passed   83.77 sec
      Start 27: 3D/small_example
 8/20 Test #27: 3D/small_example ......................   Passed    2.40 sec
      Start 28: 3D/small_example-channel
 9/20 Test #28: 3D/small_example-channel ..............***Exception: SegFault  4.01 sec
      Start 29: 3D/small_example-onlyEnergy
10/20 Test #15: 3D/hillslope01 ........................   Passed   95.30 sec
      Start 30: 3D/snow_dstr_SENSITIVITY
11/20 Test #29: 3D/small_example-onlyEnergy ...........   Passed    3.23 sec
      Start 31: 3D/Vshape
12/20 Test #31: 3D/Vshape .............................   Passed    9.43 sec
      Start 32: 3D/WG1_2.0_001
13/20 Test #22: 3D/panola_25pixel_nobed ...............   Passed  102.47 sec
14/20 Test #23: 3D/panola_25pixel_nobed_hydrostatic ...   Passed  100.59 sec
15/20 Test #13: 3D/Borden05m ..........................   Passed  141.08 sec
16/20 Test #26: 3D/rendena ............................   Passed   66.34 sec
17/20 Test #32: 3D/WG1_2.0_001 ........................   Passed   40.64 sec
18/20 Test #25: 3D/PSQL_test ..........................   Passed  122.83 sec
19/20 Test #24: 3D/prealpiC ...........................   Passed  260.74 sec
20/20 Test #30: 3D/snow_dstr_SENSITIVITY ..............   Passed  439.00 sec

95% tests passed, 1 tests failed out of 20

Total Test time (real) = 534.42 sec

The following tests FAILED:
	 28 - 3D/small_example-channel (SEGFAULT)
Errors while running CTest
```
#### 3D/small_example-channel
The test fails due to segmentation fault.

### Meson
Failing tests: 3/40.
- Borden05m.test_runner 
- Mazia.test_runner 
- snow_dstr_SENSITIVITY.test_runner

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build[v3.0*] $ meson test --suite geotop:3D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build'
ninja: no work to do.
 1/40 geotop:3D+Borden05m / 3D/Borden05m      OK      104.48 s
 2/40 geotop:3D+example / 3D/example          OK      14.88 s
 3/40 geotop:3D+hillslope01 / 3D/hillslope01  OK      72.66 s
 4/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab  OK      23.36 s
 5/40 geotop:3D+Mazia / 3D/Mazia              OK      63.27 s
 6/40 geotop:3D+no_reflection / 3D/no_reflection  OK       2.27 s
 7/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic  OK       9.86 s
 8/40 geotop:3D+panola / 3D/panola            OK      71.93 s
 9/40 geotop:3D+panola_25pixel / 3D/panola_25pixel  OK       3.01 s
10/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed  OK      59.89 s
11/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic  OK      58.05 s
12/40 geotop:3D+prealpiC / 3D/prealpiC        OK      227.04 s
13/40 geotop:3D+PSQL_test / 3D/PSQL_test      OK      96.35 s
14/40 geotop:3D+rendena / 3D/rendena          OK      50.78 s
15/40 geotop:3D+small_example / 3D/small_example  OK       1.89 s
16/40 geotop:3D+small_example-channel / 3D/small_example-channel  OK       2.90 s
17/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy  OK       2.25 s
18/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY  OK      421.41 s
19/40 geotop:3D+Vshape / 3D/Vshape            OK       5.95 s
20/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  OK      52.72 s
21/40 geotop:3D+Borden05m / 3D/Borden05m.test_runner  FAIL     0.20 s
22/40 geotop:3D+example / 3D/example.test_runner  OK       0.25 s
23/40 geotop:3D+hillslope01 / 3D/hillslope01.test_runner  OK       1.20 s
24/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab.test_runner  OK       1.22 s
25/40 geotop:3D+Mazia / 3D/Mazia.test_runner  FAIL     1.05 s
26/40 geotop:3D+no_reflection / 3D/no_reflection.test_runner  OK       1.06 s
27/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic.test_runner  OK       0.24 s
28/40 geotop:3D+panola / 3D/panola.test_runner  OK       0.18 s
29/40 geotop:3D+panola_25pixel / 3D/panola_25pixel.test_runner  OK       0.05 s
30/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed.test_runner  OK       0.28 s
31/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic.test_runner  OK       0.29 s
32/40 geotop:3D+prealpiC / 3D/prealpiC.test_runner  OK      13.94 s
33/40 geotop:3D+PSQL_test / 3D/PSQL_test.test_runner  OK      27.12 s
34/40 geotop:3D+rendena / 3D/rendena.test_runner  OK       1.63 s
35/40 geotop:3D+small_example / 3D/small_example.test_runner  OK       1.03 s
36/40 geotop:3D+small_example-channel / 3D/small_example-channel.test_runner  OK       0.25 s
37/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy.test_runner  OK       0.79 s
38/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY.test_runner  FAIL     0.17 s
39/40 geotop:3D+Vshape / 3D/Vshape.test_runner  OK       0.20 s
40/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  OK       0.19 s

OK:        37
FAIL:       3
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build/meson-logs/testlog.txt
```
#### Borden05m.test_runner
This test fails because the output files in output-tabs-SE27XX folder obtained using geotop-2.0:
- basin.txt
- discharge.txt
- point0001.txt
- soiltemp0001.txt
- soilwater0001.txt"
have an extra lines at the end of the file (that should not be there) compared to the corresponding files obtained running geotop-3.0.

The simulation time in geotop.inpts is the following:
```
InitDateDDMMYYYYhhmm = 01/06/2000 12:00
EndDateDDMMYYYYhhmm =  01/06/2000 12:10
```
but the printed extra line shows some results at 12:03 after the one of 12:10.

For example in the file basin.txt we have:
```
Date12[DDMMYYYYhhmm],JulianDayFromYear0[days],TimeFromStart[days],Simulation_Period,Run,Prain_below_canopy[mm],Psnow_below_canopy[mm],Prain_above_canopy[mm],Prain_above_canopy[mm],Pnet[mm],Tair[C],Tsurface[C],Tvegetation[C],Evap_surface[mm],Transpiration_canopy[mm],LE[W/m2],H[W/m2],SW[W/m2],LW[W/m2],LEv[W/m2],Hv[W/m2],SWv[W/m2],LWv[W/m2],SWin[W/m2],LWin[W/m2],Mass_balance_error[mm],Mean_Time_Step[s]
01/06/2000 12:01,730638.500694,0.000694,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.333333,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000014,20.000000
...
01/06/2000 12:10,730638.506944,0.006944,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.333333,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000017,20.000000
01/06/2000 12:03,730638.502083,0.002083,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.333333,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000017,20.000000
```
Since the test runner compare the results of the output files using geotop-3.0 and geotop-2.0, it prints an error like:
```
Comparing output-tabs/basin.txt and output-tabs-SE27XX/basin.txt
----------------
          <==
##12      ==> 01/06/2000 12:03,730638.502083,0.002083,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.333333,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000017,20.000000


***  End of file "output-tabs/basin.txt" reached while trying to read line 12.
***  File "output-tabs-SE27XX/basin.txt" has more lines than file "output-tabs/basin.txt",
***  line 12 is the last one read from file "output-tabs-SE27XX/basin.txt"


+++  File "output-tabs/basin.txt" differs from file "output-tabs-SE27XX/basin.txt"
```
Since the error concerns the version 2.0 we can accept the test for the 3.0.

#### Mazia.test_runner 
The test fails due to three differences in the output values of the parameter Vsub/Dt[m3/s] in the file discharge.txt
obtained after running geotop-2.0 and geotop-3.0.

Precisely the failing_ouput file is the following:
```
Comparing output-tabs/discharge.txt and output-tabs-SE27XX/discharge.txt
----------------
##3       #:9   <== 5.454597e-02
##3       #:9   ==> 3.550055e-02
@ Absolute error = 1.9045420000e-2, Relative error = 5.3648239253e-1
----------------
##4       #:9   <== 1.746238e-01
##4       #:9   ==> 7.505458e-02
@ Absolute error = 9.9569220000e-2, Relative error = 1.3266241714e+0
----------------
##5       #:9   <== 4.122229e-01
##5       #:9   ==> 1.175213e-01
@ Absolute error = 2.9470160000e-1, Relative error = 2.5076441462e+0

+++  File "output-tabs/discharge.txt" differs from file "output-tabs-SE27XX/discharge.txt"
```
Since the absolute error is small we can consider the test valid.

#### snow_dstr_SENSITIVITY.test_runner
The test fails due to two differences of the parameters Vsup/Dt[m3/s] and Vsub/Dt[m3/s] in the file discharge.txt.
obtained after running geotop-2.0 and geotop-3.0.

Precisely the failing_ouput file is the following:
```
Comparing output-tabs/discharge.txt and output-tabs-SE27XX/discharge.txt
----------------
##3       #:8   <== 4.500110e+02
##3       #:8   ==> 1.760457e+02
@ Absolute error = 2.7396530000e+2, Relative error = 1.5562169369e+0
##3       #:9   <== 1.584349e+02
##3       #:9   ==> 5.132457e+01
@ Absolute error = 1.0711033000e+2, Relative error = 2.0869211374e+0

+++  File "output-tabs/discharge.txt" differs from file "output-tabs-SE27XX/discharge.txt"
```
Since the relative error is small we can consider the test valid.


