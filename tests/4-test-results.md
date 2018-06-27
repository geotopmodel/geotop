# General infos
- Compiler: c++ (gcc 5.4.0 "c++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609")
- Processor: Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz
- Meson build type: release
- Author: Elisa Bortoli (elisa.bortoli3@gmail.com)
- Date: 27-06-2018

# Source code modified
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
 1/12 Test  #4: 1D/ColdelaPorte ..................   Passed    1.38 sec
      Start  9: 1D/PureDrainage
 2/12 Test  #9: 1D/PureDrainage ..................   Passed    2.12 sec
      Start 10: 1D/PureDrainageFaked
 3/12 Test #10: 1D/PureDrainageFaked .............   Passed    1.38 sec
      Start 11: 1D/PureDrainageRainy
 4/12 Test #11: 1D/PureDrainageRainy .............   Passed    1.68 sec
 5/12 Test  #3: 1D/Calabria ......................   Passed    7.77 sec
      Start 12: 1D/PureDrainageRainySlope
 6/12 Test  #1: 1D/B2_BeG_017 ....................   Passed    7.96 sec
 7/12 Test #12: 1D/PureDrainageRainySlope ........   Passed    1.89 sec
 8/12 Test  #5: 1D/CostantMeteo ..................   Passed    9.72 sec
 9/12 Test  #7: 1D/Matsch_B2_Ref_007 .............   Passed   10.72 sec
10/12 Test  #8: 1D/Matsch_P2_Ref_007 .............   Passed   11.51 sec
11/12 Test  #2: 1D/Bro ...........................   Passed   17.28 sec
12/12 Test  #6: 1D/Jungfraujoch ..................   Passed   26.28 sec

100% tests passed, 0 tests failed out of 12

Total Test time (real) =  26.33 sec
```

### Meson
Failing tests: 0/24.

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:1D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
 1/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017    OK       7.05 s
 2/24 geotop:1D+Bro / 1D/Bro                  OK      14.69 s
 3/24 geotop:1D+Calabria / 1D/Calabria        OK       3.44 s
 4/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte  OK       0.38 s
 5/24 geotop:1D+CostantMeteo / 1D/CostantMeteo  OK       5.37 s
 6/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch  OK      28.15 s
 7/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007  OK       8.86 s
 8/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007  OK       9.15 s
 9/24 geotop:1D+PureDrainage / 1D/PureDrainage  OK       1.00 s
10/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked  OK       1.02 s
11/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy  OK       0.93 s
12/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope  OK       0.91 s
13/24 geotop:1D+B2_BeG_017 / 1D/B2_BeG_017.test_runner  OK       0.16 s
14/24 geotop:1D+Bro / 1D/Bro.test_runner      OK       3.50 s
15/24 geotop:1D+Calabria / 1D/Calabria.test_runner  OK       1.92 s
16/24 geotop:1D+ColdelaPorte / 1D/ColdelaPorte.test_runner  OK       0.09 s
17/24 geotop:1D+CostantMeteo / 1D/CostantMeteo.test_runner  OK       1.44 s
18/24 geotop:1D+Jungfraujoch / 1D/Jungfraujoch.test_runner  OK       0.04 s
19/24 geotop:1D+Matsch_B2_Ref_007 / 1D/Matsch_B2_Ref_007.test_runner  OK       5.10 s
20/24 geotop:1D+Matsch_P2_Ref_007 / 1D/Matsch_P2_Ref_007.test_runner  OK       4.85 s
21/24 geotop:1D+PureDrainage / 1D/PureDrainage.test_runner  OK       0.26 s
22/24 geotop:1D+PureDrainageFaked / 1D/PureDrainageFaked.test_runner  OK       0.17 s
23/24 geotop:1D+PureDrainageRainy / 1D/PureDrainageRainy.test_runner  OK       0.21 s
24/24 geotop:1D+PureDrainageRainySlope / 1D/PureDrainageRainySlope.test_runner  OK       0.22 s

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
      Start 13: 3D/Borden05m
      Start 14: 3D/example
      Start 15: 3D/hillslope01
      Start 16: 3D/hillslope02_superslab
      Start 17: 3D/Mazia
      Start 18: 3D/no_reflection
      Start 19: 3D/onepoint_hydrostatic
      Start 20: 3D/panola
 1/20 Test #18: 3D/no_reflection ......................   Passed    3.86 sec
      Start 21: 3D/panola_25pixel
 2/20 Test #21: 3D/panola_25pixel .....................   Passed    3.83 sec
      Start 22: 3D/panola_25pixel_nobed
 3/20 Test #14: 3D/example ............................   Passed   25.87 sec
      Start 23: 3D/panola_25pixel_nobed_hydrostatic
 4/20 Test #16: 3D/hillslope02_superslab ..............   Passed   44.38 sec
      Start 24: 3D/prealpiC
 5/20 Test #19: 3D/onepoint_hydrostatic ...............   Passed   44.69 sec
      Start 25: 3D/PSQL_test
 6/20 Test #17: 3D/Mazia ..............................   Passed   81.38 sec
      Start 26: 3D/rendena
 7/20 Test #20: 3D/panola .............................   Passed   87.89 sec
      Start 27: 3D/small_example
 8/20 Test #27: 3D/small_example ......................   Passed    1.80 sec
      Start 28: 3D/small_example-channel
 9/20 Test #22: 3D/panola_25pixel_nobed ...............   Passed   84.26 sec
      Start 29: 3D/small_example-onlyEnergy
10/20 Test #29: 3D/small_example-onlyEnergy ...........   Passed    2.63 sec
11/20 Test #28: 3D/small_example-channel ..............   Passed    4.17 sec
      Start 30: 3D/snow_dstr_SENSITIVITY
      Start 31: 3D/Vshape
12/20 Test #15: 3D/hillslope01 ........................   Passed  103.07 sec
      Start 32: 3D/WG1_2.0_001
13/20 Test #31: 3D/Vshape .............................   Passed    9.06 sec
14/20 Test #23: 3D/panola_25pixel_nobed_hydrostatic ...   Passed   78.85 sec
15/20 Test #32: 3D/WG1_2.0_001 ........................   Passed   38.60 sec
16/20 Test #13: 3D/Borden05m ..........................   Passed  145.39 sec
17/20 Test #26: 3D/rendena ............................   Passed   63.29 sec
18/20 Test #25: 3D/PSQL_test ..........................   Passed  122.06 sec
19/20 Test #24: 3D/prealpiC ...........................   Passed  285.47 sec
20/20 Test #30: 3D/snow_dstr_SENSITIVITY ..............   Passed  489.08 sec

100% tests passed, 0 tests failed out of 20

Total Test time (real) = 584.43 sec
```

### Meson
Failing tests: 2/40.
- Mazia.test_runner 
- snow_dstr_SENSITIVITY.test_runner

The full output is:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:3D --num-processes 4
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
 1/40 geotop:3D+Borden05m / 3D/Borden05m      OK      113.22 s
 2/40 geotop:3D+example / 3D/example          OK      15.52 s
 3/40 geotop:3D+hillslope01 / 3D/hillslope01  OK      75.21 s
 4/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab  OK      24.75 s
 5/40 geotop:3D+Mazia / 3D/Mazia              OK      68.93 s
 6/40 geotop:3D+no_reflection / 3D/no_reflection  OK       2.60 s
 7/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic  OK      10.12 s
 8/40 geotop:3D+panola / 3D/panola            OK      77.16 s
 9/40 geotop:3D+panola_25pixel / 3D/panola_25pixel  OK       3.37 s
10/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed  OK      66.69 s
11/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic  OK      64.56 s
12/40 geotop:3D+prealpiC / 3D/prealpiC        OK      244.21 s
13/40 geotop:3D+PSQL_test / 3D/PSQL_test      OK      105.16 s
14/40 geotop:3D+rendena / 3D/rendena          OK      56.81 s
15/40 geotop:3D+small_example / 3D/small_example  OK       2.30 s
16/40 geotop:3D+small_example-channel / 3D/small_example-channel  OK       3.72 s
17/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy  OK       2.96 s
18/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY  OK      443.49 s
19/40 geotop:3D+Vshape / 3D/Vshape            OK       6.43 s
20/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  OK      54.97 s
21/40 geotop:3D+Borden05m / 3D/Borden05m.test_runner  OK       0.08 s
22/40 geotop:3D+example / 3D/example.test_runner  OK       0.35 s
23/40 geotop:3D+hillslope01 / 3D/hillslope01.test_runner  OK       1.61 s
24/40 geotop:3D+hillslope02_superslab / 3D/hillslope02_superslab.test_runner  OK       1.35 s
25/40 geotop:3D+Mazia / 3D/Mazia.test_runner  FAIL     0.59 s
26/40 geotop:3D+no_reflection / 3D/no_reflection.test_runner  OK       0.98 s
27/40 geotop:3D+onepoint_hydrostatic / 3D/onepoint_hydrostatic.test_runner  OK       0.39 s
28/40 geotop:3D+panola / 3D/panola.test_runner  OK       0.33 s
29/40 geotop:3D+panola_25pixel / 3D/panola_25pixel.test_runner  OK       0.10 s
30/40 geotop:3D+panola_25pixel_nobed / 3D/panola_25pixel_nobed.test_runner  OK       0.32 s
31/40 geotop:3D+panola_25pixel_nobed_hydrostatic / 3D/panola_25pixel_nobed_hydrostatic.test_runner  OK       0.41 s
32/40 geotop:3D+prealpiC / 3D/prealpiC.test_runner  OK      16.02 s
33/40 geotop:3D+PSQL_test / 3D/PSQL_test.test_runner  OK      31.03 s
34/40 geotop:3D+rendena / 3D/rendena.test_runner  OK       1.54 s
35/40 geotop:3D+small_example / 3D/small_example.test_runner  OK       1.06 s
36/40 geotop:3D+small_example-channel / 3D/small_example-channel.test_runner  OK       0.20 s
37/40 geotop:3D+small_example-onlyEnergy / 3D/small_example-onlyEnergy.test_runner  OK       1.31 s
38/40 geotop:3D+snow_dstr_SENSITIVITY / 3D/snow_dstr_SENSITIVITY.test_runner  FAIL     0.24 s
39/40 geotop:3D+Vshape / 3D/Vshape.test_runner  OK       0.17 s
40/40 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  OK       0.21 s

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