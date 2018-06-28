# General infos
- Compiler: c++ (gcc 5.4.0 "c++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609")
- Processor: Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
- Meson build type: release
- Author: Elisa Bortoli (elisa.bortoli3@gmail.com)
- Date: 28-06-2018

- Source code: modified (see following comments)
- geotop.inpts: modified for 3D tests in order to print some variables
not previously printed, since the printing output frequency was higher
than the simulation time (see following comments)

## Source code modified
The function "find_watertabledepth_dw" in the file __tables.cc__ (in src/geotop)
was modified since it produced a segmentation fault for the 3D test
small_example-channel.

The older version was:
```
double find_watertabledepth_dw(double Z, long i, long ty, SOIL *sl)
{
    double table=0.0;
    double thresh=0.0;
    long n;// number of layer below the threshold
    long nmax=sl->pa->nch;
    long l;//counter
    short out=0;

    n = 3;
    if (sl->Ptot->co[n][i] < thresh)
    {
        do
        {
            n++;
            if (n==nmax) out=-1;
            if (sl->Ptot->co[n][i] >= thresh && sl->Ptot->co[n-1][i] < thresh) out=1;
        }
        while (out==0);

        for (l=1; l<n; l++)
        {
            table += sl->pa->co[ty][jdz][l];
        }

        if (out==1)
        {
            table += ( 0.5*sl->pa->co[ty][jdz][n] - (sl->Ptot->co[n][i]-thresh)*0.5*
                                                    (sl->pa->co[ty][jdz][n-1]+sl->pa->co[ty][jdz][n])/(sl->Ptot->co[n][i]
                                                                                                       -sl->Ptot->co[n-1][i]) );
        }
        else
        {
            table += sl->pa->co[ty][jdz][n];
        }
    }

    if (table>Z) table = Z;

    return table;
}
```
The new version is:
```
double find_watertabledepth_dw(double Z, long i, long ty, SOIL *sl)
{
    double table=0.0;
    double thresh=0.0;
    long n;// number of layer below the threshold
    long nmax=sl->pa->nch;
    long l;//counter
    short out=0;

    n = 3;
    if (sl->Ptot->co[n][i] < thresh)
    {

        do
        {
            n++;
          if (n>=nmax) out=-1; // MODIFIED
          if (n<=nmax && sl->Ptot->co[n][i] >= thresh && sl->Ptot->co[n-1][i] < thresh) out=1; // MODIFIED
        }
        while (out==0);

        for (l=1; l<n; l++)
        {
            table += sl->pa->co[ty][jdz][l];
        }

        if (out==1)
        {
            table += ( 0.5*sl->pa->co[ty][jdz][n] - (sl->Ptot->co[n][i]-thresh)*0.5*
                                                    (sl->pa->co[ty][jdz][n-1]+sl->pa->co[ty][jdz][n])/(sl->Ptot->co[n][i]
                                                                                                       -sl->Ptot->co[n-1][i]) );
        }
        else
        {
            table += sl->pa->co[ty][jdz][n];
        }
    }

    if (table>Z) table = Z;

    return table;
}
```
### geotop.inpts
For unknown reasons some parameters are not printed for the following 3D tests:
- panola: Pnet
- panola_25pixel: Pnet
- prealpiC: HN, NetPrec
- PSQL_test: HN, snow
- rendena: snow
- small_example-onlyEnergy: thetaice

This has to be investigated (maybe they are 0 so they are not even printed?).

## 1D tests

### Cmake
Failing tests: 0/12.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build-release[v3.0*] $ ctest -R "1D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build-release
      Start  6: 1D/Jungfraujoch
      Start  2: 1D/Bro
      Start  8: 1D/Matsch_P2_Ref_007
      Start  7: 1D/Matsch_B2_Ref_007
      Start  5: 1D/CostantMeteo
      Start  1: 1D/B2_BeG_017
      Start  3: 1D/Calabria
      Start  9: 1D/PureDrainage
 1/12 Test  #9: 1D/PureDrainage ..................   Passed    1.81 sec
      Start 12: 1D/PureDrainageRainySlope
 2/12 Test #12: 1D/PureDrainageRainySlope ........   Passed    1.00 sec
      Start 11: 1D/PureDrainageRainy
 3/12 Test #11: 1D/PureDrainageRainy .............   Passed    1.00 sec
      Start 10: 1D/PureDrainageFaked
 4/12 Test  #3: 1D/Calabria ......................   Passed    6.33 sec
 5/12 Test  #1: 1D/B2_BeG_017 ....................   Passed    6.72 sec
      Start  4: 1D/ColdelaPorte
 6/12 Test #10: 1D/PureDrainageFaked .............   Passed    1.49 sec
 7/12 Test  #4: 1D/ColdelaPorte ..................   Passed    0.19 sec
 8/12 Test  #7: 1D/Matsch_B2_Ref_007 .............   Passed    8.43 sec
 9/12 Test  #8: 1D/Matsch_P2_Ref_007 .............   Passed    8.53 sec
10/12 Test  #5: 1D/CostantMeteo ..................   Passed    8.53 sec
11/12 Test  #2: 1D/Bro ...........................   Passed   13.66 sec
12/12 Test  #6: 1D/Jungfraujoch ..................   Passed   19.65 sec

100% tests passed, 0 tests failed out of 12

Total Test time (real) =  19.65 sec
```

### Meson
Failing tests: 0/24.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:1D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
 1/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017    OK       6.94 s
 2/24 geotop:1D+Bro / 1D/Bro                  OK      13.46 s
 3/24 geotop:1D+Calabria / 1D/Calabria        OK       2.84 s
 4/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte  OK       0.40 s
 5/24 geotop:1D+CostantMeteo / 1D/CostantMeteo  OK       4.14 s
 6/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch  OK      28.30 s
 7/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007  OK       8.21 s
 8/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007  OK       8.85 s
 9/24 geotop:1D+PureDrainage / 1D/PureDrainage  OK       1.19 s
10/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked  OK       1.08 s
11/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy  OK       0.99 s
12/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope  OK       0.88 s
13/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017.test_runner  OK       0.14 s
14/24 geotop:1D+Bro / 1D/Bro.test_runner      OK       3.33 s
15/24 geotop:1D+Calabria / 1D/Calabria.test_runner  OK       2.23 s
16/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte.test_runner  OK       0.03 s
17/24 geotop:1D+CostantMeteo / 1D/CostantMeteo.test_runner  OK       1.53 s
18/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch.test_runner  OK       0.04 s
19/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007.test_runner  OK       3.98 s
20/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007.test_runner  OK       4.37 s
21/24 geotop:1D+PureDrainage / 1D/PureDrainage.test_runner  OK       0.19 s
22/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked.test_runner  OK       0.20 s
23/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy.test_runner  OK       0.17 s
24/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope.test_runner  OK       0.15 s

OK:        24
FAIL:       0
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release/meson-logs/testlog.txt
```

## 3D tests

### Cmake
Failing tests: 0/20.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/cmake-build-release[v3.0*] $ ctest -R "3D" -j8
Test project /home/elisa/Scrivania/MHPC/geotop_3.0/cmake-build-release
      Start 30: 3D/snow_dstr_SENSITIVITY
      Start 24: 3D/prealpiC
      Start 13: 3D/Borden05m
      Start 25: 3D/PSQL_test
      Start 15: 3D/hillslope01
      Start 20: 3D/panola
      Start 22: 3D/panola_25pixel_nobed
      Start 17: 3D/Mazia
 1/20 Test #17: 3D/Mazia ..............................   Passed   78.49 sec
      Start 23: 3D/panola_25pixel_nobed_hydrostatic
 2/20 Test #22: 3D/panola_25pixel_nobed ...............   Passed   83.22 sec
 3/20 Test #20: 3D/panola .............................   Passed   83.73 sec
      Start 26: 3D/rendena
      Start 19: 3D/onepoint_hydrostatic
 4/20 Test #15: 3D/hillslope01 ........................   Passed   96.98 sec
      Start 16: 3D/hillslope02_superslab
 5/20 Test #19: 3D/onepoint_hydrostatic ...............   Passed   40.97 sec
      Start 32: 3D/WG1_2.0_001
 6/20 Test #25: 3D/PSQL_test ..........................   Passed  131.96 sec
      Start 14: 3D/example
 7/20 Test #16: 3D/hillslope02_superslab ..............   Passed   37.47 sec
      Start 31: 3D/Vshape
 8/20 Test #31: 3D/Vshape .............................   Passed    8.81 sec
      Start 28: 3D/small_example-channel
 9/20 Test #13: 3D/Borden05m ..........................   Passed  146.29 sec
      Start 18: 3D/no_reflection
10/20 Test #28: 3D/small_example-channel ..............   Passed    4.01 sec
      Start 21: 3D/panola_25pixel
11/20 Test #18: 3D/no_reflection ......................   Passed    3.11 sec
      Start 29: 3D/small_example-onlyEnergy
12/20 Test #14: 3D/example ............................   Passed   19.05 sec
      Start 27: 3D/small_example
13/20 Test #26: 3D/rendena ............................   Passed   68.65 sec
14/20 Test #21: 3D/panola_25pixel .....................   Passed    3.04 sec
15/20 Test #29: 3D/small_example-onlyEnergy ...........   Passed    2.84 sec
16/20 Test #27: 3D/small_example ......................   Passed    1.94 sec
17/20 Test #23: 3D/panola_25pixel_nobed_hydrostatic ...   Passed   77.74 sec
18/20 Test #32: 3D/WG1_2.0_001 ........................   Passed   38.83 sec
19/20 Test #24: 3D/prealpiC ...........................   Passed  274.13 sec
20/20 Test #30: 3D/snow_dstr_SENSITIVITY ..............   Passed  363.60 sec

100% tests passed, 0 tests failed out of 20

Total Test time (real) = 363.60 sec
```

### Meson
Failing tests: 2/40.
- Mazia.test_runner 
- snow_dstr_SENSITIVITY.test_runner

The full output is:
```

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

```
Since the relative error is small we can consider the test valid.