
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

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "struct.geotop.h"
//#include "input.h"
#include "constants.h"
#include "../libraries/ascii/tabs.h"
//#include "../libraries/ascii/extensions.h"
#include "../libraries/ascii/rw_maps.h"
#include "times.h"
#include "pedo.funct.h"
#include "input.h"
#include <meteoio/MeteoIO.h>
#include <string>

extern long number_novalue, number_absent;
extern char *string_novalue;

//extern char *WORKING_DIRECTORY;
extern std::string WORKING_DIRECTORY;

extern char **files;

extern long *opnt, nopnt, *obsn, nobsn, *osnw, nosnw, *oglc, noglc, *osl, nosl;
extern short *ipnt, *ibsn;
extern char **hpnt, **hbsn, **hsnw, **hglc, **hsl;
//extern char *keywords_num[num_par_number] , *keywords_char[num_par_char];
extern std::string keywords_num[], keywords_char[];


//extern char *SuccessfulRunFile, *FailedRunFile;
extern std::string SuccessfulRunFile;
extern std::string FailedRunFile;

short read_inpts_par(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, std::string filename, FILE *flog) ;

void assign_numeric_parameters(Par *par, Land *land, Times *times, Soil *sl, Meteo *met, InitTools *itools, double **num_param, long *num_param_components, std::string keyword[], FILE *flog) ;

char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, char **keyword);
char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, std::string keyword[]);

double assignation_number(FILE *f, long i, long j, char **keyword, double **num_param, long *num_param_components, double default_value, short code_error);

double assignation_number(FILE *f, long i, long j, std::string keyword[], double **num_param, long *num_param_components, double default_value, short code_error);

char *assignation_string(FILE *f, long i, char **keyword, char **string_param);
char *assignation_string(FILE *f, long i, std::string keyword[], char **string_param);

//short read_soil_parameters(char *name, char **key_header, SOIL *sl, FILE *flog);
short read_soil_parameters(char *name, InitTools *IT, Soil *sl, long bed, FILE *flog);

//short read_point_file(char *name, char **key_header, PAR *par, FILE *flog);
short read_point_file(char *name, char **key_header, Par *par, FILE *flog);

//short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name, char **key_header, FILE *flog);
short read_meteostations_file(const GeoVector<long>& i, MeteoStations *S, char *name, char **key_header, FILE *flog);

#endif
