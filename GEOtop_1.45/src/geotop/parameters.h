
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.145 'Montebello' - 8 Nov 2010
 
 Copyright (c), 2010 - Stefano Endrizzi - Geographical Institute, University of Zurich, Switzerland - stefano.endrizzi@geo.uzh.ch 
 
 This file is part of GEOtop 1.145 'Montebello'
 
 GEOtop 1.145 'Montebello' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.145 'Montebello' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */



short read_inpts_par(PAR *par, LANDCOVER *land, TIMES *times, SOIL *sl, METEO *met, INIT_TOOLS *itools, char *filename, FILE *flog);

void assign_numeric_parameters(PAR *par, LANDCOVER *land, TIMES *times, SOIL *sl, METEO *met, INIT_TOOLS *itools, double **num_param, long *num_param_components, char **keyword, FILE *flog);

char **assign_string_parameter(FILE *f, long beg, long end, char **string_param, char **keyword);

double assignation_number(FILE *f, long i, long j, char **keyword, double **num_param, long *num_param_components, double default_value, short code_error);

char *assignation_string(FILE *f, long i, char **keyword, char **string_param);

short read_soil_parameters(char *name, char **key_header, SOIL *sl, FILE *flog);

short read_point_file(char *name, char **key_header, PAR *par, FILE *flog);

short read_meteostations_file(LONGVECTOR *i, METEO_STATIONS *S, char *name, char **key_header, FILE *flog);

