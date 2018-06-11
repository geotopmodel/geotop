# General infos
Compiler: c++ (gcc 5.4.0 "c++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609")
Processor: Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
Author: Elisa Bortoli (elisa.bortoli3@gmail.com)
Date: 28-03-2018

# 1D tests

## Cmake
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/build-release-cmake[v3.0*] $ ctest -R "1D" -j4
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/build-release-cmake
      Start  1: 1D/B2_BeG_017
      Start  2: 1D/Bro
      Start  3: 1D/Calabria
      Start  4: 1D/ColdelaPorte
 1/12 Test  #4: 1D/ColdelaPorte ..................   Passed    0.29 sec
      Start  5: 1D/CostantMeteo
 2/12 Test  #3: 1D/Calabria ......................   Passed    2.72 sec
      Start  6: 1D/Jungfraujoch
 3/12 Test  #1: 1D/B2_BeG_017 ....................   Passed    4.01 sec
      Start  7: 1D/Matsch_B2_Ref_007
 4/12 Test  #5: 1D/CostantMeteo ..................   Passed    4.22 sec
      Start  8: 1D/Matsch_P2_Ref_007
 5/12 Test  #7: 1D/Matsch_B2_Ref_007 .............   Passed    5.61 sec
      Start  9: 1D/PureDrainage
 6/12 Test  #8: 1D/Matsch_P2_Ref_007 .............   Passed    6.03 sec
      Start 10: 1D/PureDrainageFaked
 7/12 Test  #9: 1D/PureDrainage ..................   Passed    0.73 sec
      Start 11: 1D/PureDrainageRainy
 8/12 Test #10: 1D/PureDrainageFaked .............   Passed    0.86 sec
      Start 12: 1D/PureDrainageRainySlope
 9/12 Test #11: 1D/PureDrainageRainy .............   Passed    1.16 sec
10/12 Test #12: 1D/PureDrainageRainySlope ........   Passed    0.70 sec
11/12 Test  #2: 1D/Bro ...........................   Passed   13.05 sec
12/12 Test  #6: 1D/Jungfraujoch ..................   Passed   17.19 sec

100% tests passed, 0 tests failed out of 12

Total Test time (real) =  19.92 sec
```

## Meson
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/build-release-meson[v3.0*] $ meson test --suite geotop:1D
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/build-release-meson'
ninja: no work to do.
 1/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017    OK      10.84 s
 2/24 geotop:1D+Bro / 1D/Bro                  OK      15.76 s
 3/24 geotop:1D+Calabria / 1D/Calabria        OK       3.76 s
 4/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte  OK       0.66 s
 5/24 geotop:1D+CostantMeteo / 1D/CostantMeteo  OK       5.79 s
 6/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch  OK      29.53 s
 7/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007  OK      10.38 s
 8/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007  OK      11.11 s
 9/24 geotop:1D+PureDrainage / 1D/PureDrainage  OK       1.50 s
10/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked  OK       1.17 s
11/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy  OK       1.49 s
12/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope  OK       1.46 s
13/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017.test_runner  OK       0.06 s
14/24 geotop:1D+Bro / 1D/Bro.test_runner      OK       3.16 s
15/24 geotop:1D+Calabria / 1D/Calabria.test_runner  OK       1.83 s
16/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte.test_runner  OK       0.01 s
17/24 geotop:1D+CostantMeteo / 1D/CostantMeteo.test_runner  OK       1.69 s
18/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch.test_runner  OK       0.02 s
19/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007.test_runner  OK       3.70 s
20/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007.test_runner  OK       4.38 s
21/24 geotop:1D+PureDrainage / 1D/PureDrainage.test_runner  OK       0.13 s
22/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked.test_runner  OK       0.13 s
23/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy.test_runner  FAIL     0.27 s
24/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope.test_runner  FAIL     0.76 s

OK:        22
FAIL:       2
SKIP:       0
TIMEOUT:    0
```

# 3D tests


## Cmake
```
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/build-release-cmake
      Start 13: 3D/Borden05m
      Start 14: 3D/example
      Start 15: 3D/hillslope01
      Start 16: 3D/hillslope02_superslab
 1/20 Test #14: 3D/example ............................   Passed   14.83 sec
      Start 17: 3D/Mazia
 2/20 Test #16: 3D/hillslope02_superslab ..............***Exception: SegFault 25.64 sec
      Start 18: 3D/no_reflection
 3/20 Test #18: 3D/no_reflection ......................***Exception: SegFault  1.71 sec
      Start 19: 3D/onepoint_hydrostatic
 4/20 Test #19: 3D/onepoint_hydrostatic ...............***Exception: SegFault  6.11 sec
      Start 20: 3D/panola
 5/20 Test #17: 3D/Mazia ..............................   Passed   54.39 sec
      Start 21: 3D/panola_25pixel
 6/20 Test #21: 3D/panola_25pixel .....................***Exception: SegFault  2.40 sec
      Start 22: 3D/panola_25pixel_nobed
 7/20 Test #15: 3D/hillslope01 ........................***Exception: SegFault 73.62 sec
      Start 23: 3D/panola_25pixel_nobed_hydrostatic
 8/20 Test #22: 3D/panola_25pixel_nobed ...............***Exception: SegFault 17.95 sec
      Start 24: 3D/prealpiC
 9/20 Test #23: 3D/panola_25pixel_nobed_hydrostatic ...***Exception: SegFault 18.04 sec
      Start 25: 3D/PSQL_test
10/20 Test #25: 3D/PSQL_test ..........................***Failed    0.46 sec
      Start 26: 3D/rendena
11/20 Test #13: 3D/Borden05m ..........................   Passed  103.56 sec
      Start 27: 3D/small_example
12/20 Test #27: 3D/small_example ......................***Exception: SegFault  1.50 sec
      Start 28: 3D/small_example-channel
13/20 Test #20: 3D/panola .............................***Exception: SegFault 71.10 sec
      Start 29: 3D/small_example-onlyEnergy
14/20 Test #29: 3D/small_example-onlyEnergy ...........***Exception: SegFault  1.70 sec
      Start 30: 3D/snow_dstr_SENSITIVITY
15/20 Test #28: 3D/small_example-channel ..............***Exception: SegFault  2.51 sec
      Start 31: 3D/Vshape
16/20 Test #31: 3D/Vshape .............................***Exception: SegFault  1.60 sec
      Start 32: 3D/WG1_2.0_001
17/20 Test #32: 3D/WG1_2.0_001 ........................   Passed   31.34 sec
18/20 Test #26: 3D/rendena ............................***Exception: SegFault 50.61 sec
19/20 Test #24: 3D/prealpiC ...........................***Exception: SegFault 78.62 sec
20/20 Test #30: 3D/snow_dstr_SENSITIVITY ..............   Passed  415.77 sec

25% tests passed, 15 tests failed out of 20

Total Test time (real) = 523.04 sec

The following tests FAILED:
	 15 - 3D/hillslope01 (SEGFAULT)
	 16 - 3D/hillslope02_superslab (SEGFAULT)
	 18 - 3D/no_reflection (SEGFAULT)
	 19 - 3D/onepoint_hydrostatic (SEGFAULT)
	 20 - 3D/panola (SEGFAULT)
	 21 - 3D/panola_25pixel (SEGFAULT)
	 22 - 3D/panola_25pixel_nobed (SEGFAULT)
	 23 - 3D/panola_25pixel_nobed_hydrostatic (SEGFAULT)
	 24 - 3D/prealpiC (SEGFAULT)
	 25 - 3D/PSQL_test (Failed)
	 26 - 3D/rendena (SEGFAULT)
	 27 - 3D/small_example (SEGFAULT)
	 28 - 3D/small_example-channel (SEGFAULT)
	 29 - 3D/small_example-onlyEnergy (SEGFAULT)
	 31 - 3D/Vshape (SEGFAULT)
Errors while running CTest
```

## Meson
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/build-release-meson[v3.0*] $ meson test --suite geotop:3D
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/build-release-meson'
ninja: no work to do.
 1/40 geotop:3D+Borden05m / 3D/Borden05m      OK      133.80 s
 2/40 geotop:3D+example / 3D/example          OK      21.01 s
 3/40 geotop:3D+hillslope01 / 3D/hillslope01  FAIL    90.78 s
 4/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab  FAIL    32.14 s
 5/40 geotop:3D+Mazia / 3D/Mazia              OK      92.31 s
 6/40 geotop:3D+no_reflection / 3D/no_reflection  FAIL     3.75 s
 7/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic  FAIL     4.70 s
 8/40 geotop:3D+panola / 3D/panola            FAIL    87.24 s
 9/40 geotop:3D+panola_25pixel / 3D/panola_25pixel  FAIL     4.44 s
10/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed  FAIL    22.81 s
11/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic  FAIL    22.95 s
12/40 geotop:3D+prealpiC / 3D/prealpiC        FAIL    109.46 s
13/40 geotop:3D+PSQL_test / 3D/PSQL_test      FAIL     0.32 s
14/40 geotop:3D+rendena / 3D/rendena          FAIL    69.09 s
15/40 geotop:3D+small_example / 3D/small_example  FAIL     2.93 s
16/40 geotop:3D+small_example-channel / 3D/small_example-channel  OK       4.57 s
17/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy  FAIL     2.57 s
18/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY  OK      435.09 s
19/40 geotop:3D+Vshape / 3D/Vshape            FAIL     1.74 s
20/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  OK      67.50 s
21/40 geotop:3D+Borden05m / 3D/Borden05m.test_runner  OK       0.02 s
22/40 geotop:3D+example / 3D/example.test_runner  OK       0.02 s
23/40 geotop:3D+hillslope01 / 3D/hillslope01.test_runner  OK       0.06 s
24/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab.test_runner  FAIL     0.07 s
25/40 geotop:3D+Mazia / 3D/Mazia.test_runner  FAIL     0.14 s
26/40 geotop:3D+no_reflection / 3D/no_reflection.test_runner  FAIL     0.02 s
27/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic.test_runner  FAIL     0.01 s
28/40 geotop:3D+panola / 3D/panola.test_runner  OK       0.02 s
29/40 geotop:3D+panola_25pixel / 3D/panola_25pixel.test_runner  OK       0.02 s
30/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed.test_runner  FAIL     0.01 s
31/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic.test_runner  FAIL     0.01 s
32/40 geotop:3D+prealpiC / 3D/prealpiC.test_runner  FAIL     0.05 s
33/40 geotop:3D+PSQL_test / 3D/PSQL_test.test_runner  OK       0.01 s
34/40 geotop:3D+rendena / 3D/rendena.test_runner  FAIL     0.01 s
35/40 geotop:3D+small_example / 3D/small_example.test_runner  FAIL     0.03 s
36/40 geotop:3D+small_example-channel / 3D/small_example-channel.test_runner  FAIL     0.16 s
37/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy.test_runner  FAIL     0.04 s
38/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY.test_runner  FAIL     0.09 s
39/40 geotop:3D+Vshape / 3D/Vshape.test_runner  FAIL     0.01 s
40/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  FAIL     0.05 s

OK:        12
FAIL:      28
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/build-release-meson/meson-logs/testlog.txt
```
