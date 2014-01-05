
/* STATEMENT:
 overloaded function
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
    
#ifndef OUTPUT_H
#define OUTPUT_H

#include "struct.geotop.h"
#include "pedo.funct.h"
#include "../libraries/ascii/rw_maps.h"
#include "constants.h"
#include "times.h"
#include "energy.balance.h"
#include "input.h"
#include "../libraries/ascii/tabs.h"
#include "vegetation.h"
#include "tables.h"
#include "snow.h"


#include <time.h>

//void write_output(TIMES *times,WATER *wat,CHANNEL *cnet,PAR *par,TOPO *top,LAND *land,SOIL *sl,Energy *egy,SNOW *snow,GLACIER *glac,METEO *met);
  void write_output(Times *times,Water *wat,Channel *cnet,Par *par,Topo *top,Land *land,Soil *sl,Energy *egy,Snow *snow,Glacier *glac,Meteo *met);

//void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, Energy *egy, SNOW *snow, GLACIER *glac);

//void write_output_headers(long n, TIMES *times, WATER *wat, PAR *par, TOPO *top, LAND *land, SOIL *sl, Energy *egy, Snow *snow, Glacier *glac);
  void write_output_headers(long n, Times *times, Water *wat, Par *par, Topo *top, Land *land, Soil *sl, Energy *egy, Snow *snow, Glacier *glac);

//void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, DOUBLEVECTOR *n, SOIL *sl, PAR *par, double psimin, double cosslope);
  void write_soil_output(long i, long iname, double init_date, double JDfrom0, double JD, long day, long month, long year, long hour, long minute, const GeoVector<double>& n, Soil *sl, Par *par, double psimin, double cosslope);

//void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, double *var, DOUBLEVECTOR *n, double *dz, double cosslope);
  void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, double *var, const GeoVector<double>& n, double *dz, double cosslope);
  void write_soil_file(long lmin, long i, FILE *f, long d, long m, long y, long h, long mi, double JDfrom0, double JDfrom0init, const GeoMatrix<double>& var, long row, const GeoVector<double>& n, const GeoTensor<double>& dz, double cosslope);

//void write_soil_header(FILE *f, DOUBLEVECTOR *n, double *dz);
  void write_soil_header(FILE *f, const GeoVector<double>& n, const GeoTensor<double>& dz);

//void plot(char *name, long i_plot, DOUBLEVECTOR *V, short format, long **J);
void plot(std::string name, long i_plot, const GeoVector<double>& V, short format, long **J);

//double interpolate_soil(long lmin, double h, long max, double *Dz, double *Q);
  double interpolate_soil(long lmin, double h, long max, const GeoTensor<double>& Dz, const GeoMatrix<double>& Q, const long& row);

//double interpolate_soil2(long lmin, double h, long max, double *Dz, DOUBLEMATRIX *Q, long i);
  double interpolate_soil2(long lmin, double h, long max, const GeoTensor<double>& Dz, GeoMatrix<double>& Q, long i);

void write_tensorseries_soil(long lmin, std::string suf, std::string filename, short type, short format, GeoMatrix<double>& T, const GeoVector<double>& n, long **J, GeoMatrix<long>& RC, GeoTensor<double>& dz, GeoMatrix<double>& slope, short vertical);
//void fill_output_vectors(double Dt, double W, Energy *egy, SNOW *snow, GLACIER *glac, WATER *wat, METEO *met, PAR *par, TIMES *time, TOPO *top);
  void fill_output_vectors(double Dt, double W, Energy *egy, Snow *snow, Glacier *glac, Water *wat, Meteo *met, Par *par, Times *time, Topo *top);

#endif
