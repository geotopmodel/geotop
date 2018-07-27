# Found bug in geotop-2.0 version!
Because of some problems running the version 3.0 in DEBUG mode for a 3D test, we found a bug of version 2.0.
Read it!

RELEASE build type:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-release[v3.0*] $ meson test --suite geotop:WG1_2.0_001
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release'
ninja: no work to do.
1/2 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  OK      48.62 s
2/2 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  OK       0.12 s

OK:         2
FAIL:       0
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-release/meson-logs/testlog.txt
```

Same test but in DEBUG build type:
```
elisa@elisa-N552VW ~/Scrivania/MHPC/geotop_3.0/meson-build-debug[v3.0*] $ meson test --suite geotop:WG1_2.0_001
ninja: Entering directory `/home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-debug'
ninja: no work to do.
1/2 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001  FAIL    114.22 s
2/2 geotop:3D+WG1_2.0_001 / 3D/WG1_2.0_001.test_runner  FAIL     0.13 s

OK:         0
FAIL:       2
SKIP:       0
TIMEOUT:    0

Full log written to /home/elisa/Scrivania/MHPC/geotop_3.0/meson-build-debug/meson-logs/testlog.txt
```

## Problem
In __constants.h__ the parameter ```z_transp``` was defined as follows:
```
#define z_transp 10000 //soil depth responsable for canopy transpiration [mm]
```

In __input.cc__ the parameters ```z``` (total soil depth) and ```l``` (index)
were calculated as follows:
```
// Check that there aren't cell with an undefined land use value
    z = 0.;
    l = 0;

    do
    {
        l++;
        z += sl->pa->co[1][jdz][l];
    }
    while (l<Nl && z < z_transp);
```
and then the matrix root_fraction was initialized as:
```
land->root_fraction.reset(new Matrix<double>{par->n_landuses, l});
```

Then in __output.cc__ we wanted to print the root fraction of each layer as follows:
```
                for (l=1; l<=Nl; l++)
                {
                    fprintf(f," The root fraction [-] of layer %ld: %f\n",l, (*land->root_fraction)(lu,l));
                }
```

The problem occured when __z > z_transp__ since at a certain point we wanted to access the root_fraction matrix
using an l index > its column dimension.

More precisely, the involved test cases are:
- RottCatchment
- UpperAmmerCatchment
- WG1_2.0_001

### Example: WG1_2.0_001
In this example we have:
- 10 landcover type (```par->n_landuses```)
- 9 layers (```Nl```).

The total soil depth is the sum of the values in the geotop.inpts, leading to 47030 mm:
```
SoilLayerThicknesses			=	30,100,200,500,1200,3000,5000,12000,25000
```
BUT ```z_transp = 10000```, located in the layer 7, so at a certain point, accessing:
```
(*land->root_fraction)(lu,8)
```
we got an exception since the matrix ```land->root_fraction``` had only 7 columns.


## Solution
Different solutions are possible:
1. change the end index of the for loop from ```Nl``` to a ```lmax = land->root_fraction->n_col```

```
                int lmax = land->root_fraction->n_col;
                for (l=1; l<=land->root_fraction->n_col; l++)
                {
                    fprintf(f," The root fraction [-] of layer %ld: %f\n",l, (*land->root_fraction)(lu,l));
                }
```

2. increase the ```z_transp``` parameter in constants.h file
```
+ quick solution
- problems if we insert a test cases having again z > z_transp
```
3. insert a check like the following:
```
 if (z > z_transp){
        l = Nl;
        z_transp = z;
    }
```
```
+ more general solution
- a bit dirty way
```

4. modify the function root defined in __vegetation.cc__
```
+ general solution
- need time to discuss
```
The __method 1__ will be implemented.

### Optional: implementation of solution 2
The parameters z_transp was changed from a macro-defined to a global variable defined in input.cc
to be able to modify its value.
No more segmentation fault occured BUT of course the comparison with the 2.0 version failed.
I think that the version 2.0 should be corrected and all the tests should be re-runned.

Actually looking at some output file using the 2 versions, it's clear that the 2.0 is wrong
since the sum of the root fraction of all layers should be 1 (100%) but this doesn't happen.

#### RottCatchment (pointall_info_0001.txt)
```
 Land use number is 1 
 The root fraction [-] of layer 1: 0.181818
 The root fraction [-] of layer 2: 0.181818
 The root fraction [-] of layer 3: 0.181818
 The root fraction [-] of layer 4: 0.181818
 The root fraction [-] of layer 5: 0.272727
 The root fraction [-] of layer 6: 0.000000
 The root fraction [-] of layer 7: 0.000000
 The root fraction [-] of layer 8: 0.000000
 The root fraction [-] of layer 9: 0.000000
 The root fraction [-] of layer 10: 0.000000
 The root fraction [-] of layer 11: 0.000000
 The root fraction [-] of layer 12: 0.000000
 The root fraction [-] of layer 13: 0.076923
```
#### UpperAmmerCatchment (pointall_info_0001.txt)
```
 Land use number is 12 
 The root fraction [-] of layer 1: 0.222222
 The root fraction [-] of layer 2: 0.222222
 The root fraction [-] of layer 3: 0.222222
 The root fraction [-] of layer 4: 0.222222
 The root fraction [-] of layer 5: 0.111111
 The root fraction [-] of layer 6: 0.000000
 The root fraction [-] of layer 7: 0.000000
 The root fraction [-] of layer 8: 0.000000
 The root fraction [-] of layer 9: 0.000000
 The root fraction [-] of layer 10: 0.000000
 The root fraction [-] of layer 11: 0.000000
 The root fraction [-] of layer 12: 0.000000
 The root fraction [-] of layer 13: 1.000000
 ```
## Talking with Giacomo Bertoldi
We thought to have a triangular distribution of root across depth, with greater values for the first layers (the shallowest ones)
but if we check the previous results we can see that it's not like this.

### root function
The parameter ```root_fraction``` was then calculated in __vegetation.cc__
as follows:
```
void root(long n, double d, double slope, double *D, MatrixRow<double> &&root_fraction)
{

  //n = number of soil layers (from the surface) affected by root absorption
  //d = root depth (vertical) [mm]
  //slope = max lateral slope [rad]
  //D[] = soil layer thickness [mm]
  //root_fraction[] = weighting factor assigned to each layer for transpiration, sum of all root_fraction components gives 1

  long l;
  double z=0.0;
  double d_corr=d*cos(slope); //slope depth taken as ortogonal to layer boundary

  for (l=1; l<=n; l++)
    {
      z += D[l];
      if ( d_corr > z )
        {
          root_fraction[l] = D[l]/d_corr;
        }
      else
        {
          if ( d_corr > z-D[l] )
            {
              root_fraction[l] = ( d_corr - (z-D[l]) ) / d_corr;
            }
          else
            {
              root_fraction[l] = 0.0;
            }
        }
    }

}
```



