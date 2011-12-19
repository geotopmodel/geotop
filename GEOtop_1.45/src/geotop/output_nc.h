
/* STATEMENT:
 
 
 
 */
    
#ifdef USE_NETCDF

void write_output_nc(ALLDATA*);

int ncgt_var_update(void *m, void * m0, double Dt, short nlimdim, double novalue);

int ncgt_var_set_to_zero(void * m0, short nlimdim, double novalue);

#endif
