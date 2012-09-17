
/* STATEMENT:
 
 
 
 */
    
#ifdef USE_NETCDF

#include <sys/stat.h>
#include "struct.geotop.h"
#ifdef USE_HPC
#include "hpc.geotop.h"
#endif
#include "input.h"
#include "output.h"
#include "times.h"
#include "constants.h"
#include "energy.balance.h"
#include "meteo.h"
#include "water.balance.h"
#include "snow.h"
#include "blowingsnow.h"
#include "../libraries/ascii/tabs.h"
#include "deallocate.h"

#include "../gt_utilities/gt_utilities.h"
#include "../gt_utilities/gt_symbols.h"
#include "../gt_utilities/ncgt_output.h"

#include <time.h>

extern long i_sim,i_run;
extern long number_novalue;

void write_output_nc(ALLDATA*);

int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);

void deallocate_output_nc(OUTPUT_NCDATA* outnc);

void set_output_nc(ALLDATA *all);

#endif
