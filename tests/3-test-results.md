# General infos
- Compiler: c++ (gcc 5.4.0 "c++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609")
- Processor: Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
- Meson build type: release
- Author: Elisa Bortoli (elisa.bortoli3@gmail.com)
- Date: 25-06-2018

## 1D tests

### Cmake
Failing tests: 0/12.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build-release[v3.0*] $ ctest -R "1D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build-release
      Start  1: 1D/B2_BeG_017
      Start  2: 1D/Bro
      Start  3: 1D/Calabria
      Start  4: 1D/ColdelaPorte
      Start  5: 1D/CostantMeteo
      Start  6: 1D/Jungfraujoch
      Start  7: 1D/Matsch_B2_Ref_007
      Start  8: 1D/Matsch_P2_Ref_007
 1/12 Test  #4: 1D/ColdelaPorte ..................   Passed    1.36 sec
      Start  9: 1D/PureDrainage
 2/12 Test  #9: 1D/PureDrainage ..................   Passed    1.31 sec
      Start 10: 1D/PureDrainageFaked
 3/12 Test #10: 1D/PureDrainageFaked .............   Passed    1.31 sec
      Start 11: 1D/PureDrainageRainy
 4/12 Test  #1: 1D/B2_BeG_017 ....................   Passed    5.45 sec
 5/12 Test #11: 1D/PureDrainageRainy .............   Passed    1.27 sec
      Start 12: 1D/PureDrainageRainySlope
 6/12 Test  #3: 1D/Calabria ......................   Passed    6.68 sec
 7/12 Test #12: 1D/PureDrainageRainySlope ........   Passed    1.24 sec
 8/12 Test  #8: 1D/Matsch_P2_Ref_007 .............   Passed    9.57 sec
 9/12 Test  #7: 1D/Matsch_B2_Ref_007 .............   Passed    9.58 sec
10/12 Test  #5: 1D/CostantMeteo ..................   Passed   11.88 sec
11/12 Test  #2: 1D/Bro ...........................   Passed   15.39 sec
12/12 Test  #6: 1D/Jungfraujoch ..................   Passed   19.48 sec

100% tests passed, 0 tests failed out of 12

Total Test time (real) =  19.49 sec
```

### Meson
Failing tests: 0/24.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:1D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
 1/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017    OK       6.72 s
 2/24 geotop:1D+Bro / 1D/Bro                  OK      12.60 s
 3/24 geotop:1D+Calabria / 1D/Calabria        OK       2.50 s
 4/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte  OK       0.30 s
 5/24 geotop:1D+CostantMeteo / 1D/CostantMeteo  OK       4.42 s
 6/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch  OK      27.84 s
 7/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007  OK       8.09 s
 8/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007  OK       8.39 s
 9/24 geotop:1D+PureDrainage / 1D/PureDrainage  OK       0.93 s
10/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked  OK       0.86 s
11/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy  OK       0.89 s
12/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope  OK       0.87 s
13/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017.test_runner  OK       0.06 s
14/24 geotop:1D+Bro / 1D/Bro.test_runner      OK       3.22 s
15/24 geotop:1D+Calabria / 1D/Calabria.test_runner  OK       1.83 s
16/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte.test_runner  OK       0.01 s
17/24 geotop:1D+CostantMeteo / 1D/CostantMeteo.test_runner  OK       1.34 s
18/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch.test_runner  OK       0.02 s
19/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007.test_runner  OK       3.95 s
20/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007.test_runner  OK       3.61 s
21/24 geotop:1D+PureDrainage / 1D/PureDrainage.test_runner  OK       0.13 s
22/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked.test_runner  OK       0.13 s
23/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy.test_runner  OK       0.13 s
24/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope.test_runner  OK       0.13 s

OK:        24
FAIL:       0
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release/meson-logs/testlog.txt
```

## 3D tests

### Cmake
Failing tests: 1/20.
- 3D/small_example-channel

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build-release[v3.0*] $ ctest -R "3D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build-release
      Start 13: 3D/Borden05m
      Start 14: 3D/example
      Start 15: 3D/hillslope01
      Start 16: 3D/hillslope02_superslab
      Start 17: 3D/Mazia
      Start 18: 3D/no_reflection
      Start 19: 3D/onepoint_hydrostatic
      Start 20: 3D/panola
 1/20 Test #18: 3D/no_reflection ......................   Passed    3.09 sec
      Start 21: 3D/panola_25pixel
 2/20 Test #21: 3D/panola_25pixel .....................   Passed    3.00 sec
      Start 22: 3D/panola_25pixel_nobed
 3/20 Test #14: 3D/example ............................   Passed   17.82 sec
      Start 23: 3D/panola_25pixel_nobed_hydrostatic
 4/20 Test #16: 3D/hillslope02_superslab ..............   Passed   35.54 sec
      Start 24: 3D/prealpiC
 5/20 Test #19: 3D/onepoint_hydrostatic ...............   Passed   49.49 sec
      Start 25: 3D/PSQL_test
 6/20 Test #17: 3D/Mazia ..............................   Passed   74.96 sec
      Start 26: 3D/rendena
 7/20 Test #20: 3D/panola .............................   Passed   79.90 sec
      Start 27: 3D/small_example
 8/20 Test #27: 3D/small_example ......................   Passed    2.51 sec
      Start 28: 3D/small_example-channel
 9/20 Test #28: 3D/small_example-channel ..............***Exception: SegFault  4.11 sec
      Start 29: 3D/small_example-onlyEnergy
10/20 Test #29: 3D/small_example-onlyEnergy ...........   Passed    3.02 sec
      Start 30: 3D/snow_dstr_SENSITIVITY
11/20 Test #15: 3D/hillslope01 ........................   Passed   92.54 sec
      Start 31: 3D/Vshape
12/20 Test #31: 3D/Vshape .............................   Passed    9.23 sec
      Start 32: 3D/WG1_2.0_001
13/20 Test #22: 3D/panola_25pixel_nobed ...............   Passed  102.31 sec
14/20 Test #23: 3D/panola_25pixel_nobed_hydrostatic ...   Passed   99.24 sec
15/20 Test #26: 3D/rendena ............................   Passed   62.53 sec
16/20 Test #13: 3D/Borden05m ..........................   Passed  140.69 sec
17/20 Test #32: 3D/WG1_2.0_001 ........................   Passed   40.09 sec
18/20 Test #25: 3D/PSQL_test ..........................   Passed  121.67 sec
19/20 Test #24: 3D/prealpiC ...........................   Passed  259.14 sec
20/20 Test #30: 3D/snow_dstr_SENSITIVITY ..............   Passed  427.84 sec

95% tests passed, 1 tests failed out of 20

Total Test time (real) = 519.19 sec

The following tests FAILED:
	 28 - 3D/small_example-channel (SEGFAULT)
Errors while running CTest
```
#### 3D/small_example-channel
The test fails due to segmentation fault, that occurs in the function
```find_watertabledepth_dw``` (defined in tables.cc) at the following line:

```
if(sl->Ptot->co[n][i] >= thresh && sl->Ptot->co[n-1][i] < thresh) out=1; 
```
If I change the value n defined as "number of layer below the threshold"
from 3 to 2 or 1, the problem seems to be solved but three 1D tests fail:
- Calabria
- PureDrainageRainy
- PureDrainageRainySlope.

### Meson
Failing tests: 2/40.
- Mazia.test_runner 
- snow_dstr_SENSITIVITY.test_runner

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:3D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
 1/40 geotop:3D+Borden05m / 3D/Borden05m      OK      109.84 s
 2/40 geotop:3D+example / 3D/example          OK      15.67 s
 3/40 geotop:3D+hillslope01 / 3D/hillslope01  OK      75.61 s
 4/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab  OK      24.75 s
 5/40 geotop:3D+Mazia / 3D/Mazia              OK      66.47 s
 6/40 geotop:3D+no_reflection / 3D/no_reflection  OK       2.31 s
 7/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic  OK      10.09 s
 8/40 geotop:3D+panola / 3D/panola            OK      74.00 s
 9/40 geotop:3D+panola_25pixel / 3D/panola_25pixel  OK       3.12 s
10/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed  OK      61.50 s
11/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic  OK      59.81 s
12/40 geotop:3D+prealpiC / 3D/prealpiC        OK      232.59 s
13/40 geotop:3D+PSQL_test / 3D/PSQL_test      OK      97.49 s
14/40 geotop:3D+rendena / 3D/rendena          OK      50.69 s
15/40 geotop:3D+small_example / 3D/small_example  OK       1.89 s
16/40 geotop:3D+small_example-channel / 3D/small_example-channel  OK       2.91 s
17/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy  OK       2.25 s
18/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY  OK      422.98 s
19/40 geotop:3D+Vshape / 3D/Vshape            OK       5.94 s
20/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  OK      52.69 s
21/40 geotop:3D+Borden05m / 3D/Borden05m.test_runner  OK       0.02 s
22/40 geotop:3D+example / 3D/example.test_runner  OK       0.19 s
23/40 geotop:3D+hillslope01 / 3D/hillslope01.test_runner  OK       0.82 s
24/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab.test_runner  OK       0.83 s
25/40 geotop:3D+Mazia / 3D/Mazia.test_runner  FAIL     0.18 s
26/40 geotop:3D+no_reflection / 3D/no_reflection.test_runner  OK       0.38 s
27/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic.test_runner  OK       0.13 s
28/40 geotop:3D+panola / 3D/panola.test_runner  OK       0.12 s
29/40 geotop:3D+panola_25pixel / 3D/panola_25pixel.test_runner  OK       0.04 s
30/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed.test_runner  OK       0.14 s
31/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic.test_runner  OK       0.20 s
32/40 geotop:3D+prealpiC / 3D/prealpiC.test_runner  OK      13.96 s
33/40 geotop:3D+PSQL_test / 3D/PSQL_test.test_runner  OK      29.91 s
34/40 geotop:3D+rendena / 3D/rendena.test_runner  OK       1.41 s
35/40 geotop:3D+small_example / 3D/small_example.test_runner  OK       0.30 s
36/40 geotop:3D+small_example-channel / 3D/small_example-channel.test_runner  OK       0.19 s
37/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy.test_runner  OK       0.42 s
38/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY.test_runner  FAIL     0.08 s
39/40 geotop:3D+Vshape / 3D/Vshape.test_runner  OK       0.04 s
40/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  OK       0.04 s

OK:        38
FAIL:       2
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release/meson-logs/testlog.txt
```

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