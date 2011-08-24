
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

    
typedef struct {
	double snow0;
	double Tsnow0;
	double agesnow0;
    double rhosnow0;	
	double rhoglac0;
	double Dglac0;
	double Tglac0;
	char **met_col_names;
	char **soil_col_names;
	char **horizon_col_names;
	char **point_col_names;
	char **lapserates_col_names;
	char **meteostations_col_names;
	DOUBLEMATRIX *LU;
} INIT_TOOLS;



void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, 
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times);

void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par, FILE *flog);

void read_optionsfile_point(PAR *par, TOPO *top, LAND *land, SOIL *sl, TIMES *times, INIT_TOOLS *IT, FILE *flog);

void set_bedrock(SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top, DOUBLEMATRIX *LC, FILE *flog);

DOUBLETENSOR *find_Z_of_any_layer(DOUBLEMATRIX *Zsurface, DOUBLEMATRIX *slope, DOUBLEMATRIX *LC, SOIL *sl);

short file_exists(short key, FILE *flog);

double peat_thickness(double dist_from_channel);
