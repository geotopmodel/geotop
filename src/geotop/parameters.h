
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 2.0.0

 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */




short read_inpts_par(PAR *par, LAND *land, TIMES *times, SOIL *sl, METEO *met,
                     INIT_TOOLS *itools, char *filename, FILE *flog);

void assign_numeric_parameters(PAR *par, LAND *land, TIMES *times, SOIL *sl,
                               METEO *met, INIT_TOOLS *itools, double **num_param,
                               long *num_param_components, char **keyword, FILE *flog);

char **assign_string_parameter(FILE *f, long beg, long end,
                               char **string_param, char **keyword);

double assignation_number(FILE *f, long i, long j, char **keyword,
                          double **num_param, long *num_param_components, double default_value,
                          short code_error);

char *assignation_string(FILE *f, long i, char **keyword,
                         char **string_param);

short read_soil_parameters(char *name, INIT_TOOLS *IT, SOIL *sl, long bed,
                           FILE *flog);

short read_point_file(char *name, char **key_header, PAR *par, FILE *flog);

short read_meteostations_file(Vector<long> *i, METEO_STATIONS *S, char *name,
                              char **key_header, FILE *flog);

