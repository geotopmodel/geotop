
/* STATEMENT:
 
 
 
 */
#ifndef OUTPUT_NC_H
#define OUTPUT_NC_H
#ifdef USE_NETCDF
#include <sys/stat.h>
#include "struct.geotop.h"
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

//#include "../gt_utilities/gt_utilities.h"
//#include "../gt_utilities/gt_symbols.h"
//#include "../gt_utilities/ncgt_output.h"
#include "../netCDF/netcdfIO.h"
#include "../netCDF/gt_symbols.h"
#include <time.h>

extern long i_sim,i_run;
extern long number_novalue;


void write_output_nc(AllData*);

int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);

void deallocate_output_nc(OutputNCData* outnc);

void set_output_nc(AllData *all);

#endif
#endif
