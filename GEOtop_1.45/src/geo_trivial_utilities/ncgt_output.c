#ifdef USE_NETCF_ONGOING

int add_2Dmap(int ncid, DOUBLEMATRIX *m, double time, const char* dimension_time, long *k, short reinitialize, double number_novalue){
	/* define the temporal counter*/

	long a=*k;
	ncgt_put_double_vs_time(time, dimension_time, a, ncid, dimension_time);

	/* rotate map and put to netCDF */
	ncgt_put_rotate180_y_doublematrix_vs_time(m, &k, ncid, dimension_time, NC_GEOTOP_XLAT, NC_GEOTOP_YLON);
	if(reinitialize==1) initmatrix(0.0, m, m, number_novalue);
	/* 26.11.20011 to do list:
	 * put the attributes
	 *
	 *  */
	a++;
	*k=a;
}
#endif
