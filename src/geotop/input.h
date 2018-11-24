#ifndef _GEOTOP_INPUT_H
#define _GEOTOP_INPUT_H


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


struct INIT_TOOLS
{
    double swe0;
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
    std::unique_ptr<Matrix<double>> bed;
    std::unique_ptr<Tensor<double>> pa_bed;
    std::unique_ptr<Vector<double>> init_water_table_depth;
};



void get_all_input(long argc, char *argv[], TOPO *top, SOIL *sl, LAND *land,
                   METEO *met, WATER *wat, CHANNEL *cnet,
                   PAR *par, ENERGY *egy, SNOW *snow, GLACIER *glac, TIMES *times);

void read_inputmaps(TOPO *top, LAND *land, SOIL *sl, PAR *par, INIT_TOOLS *IT);

void read_optionsfile_point(PAR *par, TOPO *top, LAND *land, SOIL *sl, TIMES *times, INIT_TOOLS *IT);

void set_bedrock(INIT_TOOLS *IT, SOIL *sl, CHANNEL *cnet, PAR *par, TOPO *top, Matrix<double> *LC);

std::unique_ptr<Tensor<double>> find_Z_of_any_layer(Matrix<double> *Zsurface, Matrix<double> *slope,
                                                    Matrix<double> *LC, SOIL *sl, short point);

short file_exists(short key);

void copy_soil_state(SOIL_STATE *from, SOIL_STATE *to);

void initialize_veg_state(STATE_VEG *V, long n);

void copy_veg_state(STATE_VEG *from, STATE_VEG *to);



#endif
