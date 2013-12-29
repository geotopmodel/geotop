
/* STATEMENT:
 
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
#ifndef DEALLOCATE_H
#define DEALLOCATE_H
#include "struct.geotop.h"
#include "constants.h"
#include "snow.h"
//#include "../libraries/ascii/init.h"

extern std::vector<std::string> files;
extern FILE *ffbas, *ffpoint, *ffT, *ffTav, *ffpsi, *ffpsitot, *ffliq, *ffliqav, *ffice, *fficeav, *ffsnow, *ffglac;
extern double **odpnt, **odp, *odbsn, *odb;
extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw;
extern long *oglc, noglc, *osl, nosl;
extern short *ipnt, *ibsn;
extern char **hpnt, **hbsn, **hsnw, **hglc, **hsl;
//extern char *WORKING_DIRECTORY;
extern std::string WORKING_DIRECTORY;

extern char *string_novalue;
extern long number_novalue;
extern long Nl, Nr, Nc;
//extern T_INIT *UV;
extern TInit *UV;


//void dealloc_all(TOPO *top,SOIL *sl,LAND *land,WATER *wat,CHANNEL *cnet,PAR *par,Energy *egy,SNOW *snow, GLACIER *glac, METEO *met, TIMES *times);
  void dealloc_all(Topo *top,Soil *sl,Land *land,Water *wat,Channel *cnet,Par *par,Energy *egy,Snow *snow, Glacier *glac, Meteo *met, Times *times);
//void dealloc_meteostations(METEO_STATIONS *st);
  void dealloc_meteostations(MeteoStations *st);
//void deallocate_soil_state(SOIL_STATE *S);
  void deallocate_soil_state(SoilState *S);
//void deallocate_veg_state(STATE_VEG *V);
  void deallocate_veg_state(StateVeg *V);
//void reset_to_zero(PAR *par, SOIL *sl, LAND *land, SNOW *snow, GLACIER *glac, Energy *egy, METEO *met, WATER *wat);
  void reset_to_zero(Par *par, Soil *sl, Land *land, Snow *snow, Glacier *glac, Energy *egy, Meteo *met, Water *wat);

#endif
