
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
    
typedef struct {
	DOUBLEMATRIX *land_classes;
	DOUBLEMATRIX *met;
	double u0;
	double D;
	double snow0;
	double Tsnow0;
	double agesnow0;
    double rhosnow0;	
	double rhoglac0;
	double Dglac0;
	double Tglac0;
	STRINGBIN *met_col_names;
	DOUBLEMATRIX *LU;
} INIT_TOOLS;



/****************************************************************************************************/
/* Subroutine which get all the input file and put the variables in the apposite structs            */
/****************************************************************************************************/
void get_all_input(int argc, char *argv[], TOPO *top, SOIL *sl, LAND *land, METEO *met, WATER *wat, CHANNEL *cnet, 
					PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times);




/****************************************************************************************************/
/* Subroutine to find the the coefficents to diffuse the channel-flows ("fraction_spread" matrix):  */
/****************************************************************************************************/
DOUBLEMATRIX * De_Saint_Venant(DOUBLEVECTOR *s0,double u0,double D,double Dt);




/****************************************************************************************************/
/* year_avg_temp: calcola la temprarura media e l'escursione annuale dell'aria                   */
/* Input:	data_meteo: matrice dati meteorolgici                                                   */
/*  		intervallo_escursione [gg]: intervallo in cui viene mediato il dato di temperatura per  */
/*	    	calcolare l'escursione media annuale                                                    */
/* Output:	temp_med_ann: temperatura media annuale                                                 */
/*		    delta_temp_ann: escursione annuale di temperatura                                       */
/****************************************************************************************************/
void year_avg_temp(DOUBLEVECTOR *Tdata, double ndays, double *T, double *DT, double Dt);

long row(double N, long nrows, long i, T_INIT *UV);

long col(double E, long ncols, long i, T_INIT *UV);

void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par);

void read_parameterfile(char *name, PAR *par, INIT_TOOLS *itools);

void read_optionsfile_distr(char *name, PAR *par, TIMES *times, DOUBLEMATRIX *Z0);

void read_optionsfile_point(char *name, PAR *par, TOPO *top, LAND *land, SOIL *sl, INIT_TOOLS *IT);

double **read_horizon(char *name, long i);

void init_meteo_stations(DOUBLEMATRIX *INPUTmeteo, METEO_STATIONS *st);

void ReadMeteoHeader(FILE *f, STRINGBIN *ColDescr, long offset, long *ncols, long *MeteoCont);

double **read_datameteo(FILE *f, long offset, long ncols, double ndef);

void read_inpts_par(PAR *par, TIMES *times, char *program, char *ext);

DOUBLETENSOR *read_soil_parameters(char *name);

DOUBLEMATRIX *depitted(SHORTMATRIX *DD, DOUBLEMATRIX *Z);

void assign_recovered(char *name, double **assign, PAR *par, DOUBLEMATRIX *Z1, DOUBLEMATRIX *Z2);

void assign_recovered_long(char *name, long **assign, PAR *par, DOUBLEMATRIX *Z1, DOUBLEMATRIX *Z2);

long n_entries_matrix(TOPO *top, DOUBLEMATRIX *LC, long n);

//void complete_Ap_Ai(TOPO *top, DOUBLEMATRIX *LC, long n);

void find_max_constraint( DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long *R, long *C);

void next_down_channel_pixel( long r, long c, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH, long *R, long *C);

short neighboring_down_channel_pixel( long r, long c, long ir, long ic, DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, SHORTMATRIX *pixel_type, LONGMATRIX *CH);

void i_lrc_cont(DOUBLEMATRIX *LC, long ***i, LONGMATRIX *lrc);

void j_rc_cont(DOUBLEMATRIX *LC, long **j, LONGMATRIX *rc);

void searchDD(long **J, LONGMATRIX *RC, SHORTMATRIX *DD, DOUBLEMATRIX *LC, LONGMATRIX *up, LONGVECTOR *down);

//void cont_nonzero_values_matrix(WATER *wat, DOUBLEMATRIX *LC, long ***i, long n, short point);

void cont_nonzero_values_matrix2(long *tot, DOUBLEMATRIX *LC, LONGMATRIX *lrc, long ***i, long n, short point);

void cont_nonzero_values_matrix3(LONGVECTOR *Lp, LONGVECTOR *Li, DOUBLEMATRIX *LC, LONGMATRIX *lrc, long ***i, long n, short point);

void cont_nonzero_values_matrix4(LONGVECTOR *Lp, LONGVECTOR *Li, LONGVECTOR *Up, LONGVECTOR *Ui, DOUBLEMATRIX *LC, LONGMATRIX *lrc, long ***i, long n, short point);

void remove_outdraining_subbasins(SHORTMATRIX *DD, SHORTMATRIX *net, DOUBLEMATRIX *LC, double novalue);

double peat_thickness(double dist_from_channel);

