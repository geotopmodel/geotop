/* STATEMENT:

GEO_TOP MODELS THE ENERGY AND WATER FLUXES AT LAND SURFACE
GEOtop-Version 0.9375-Subversion Mackenzie

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon, Emanuele Cordano, Matteo Dall'Amico

 LICENSE:

 This file is part of GEOtop 0.9375 Mackenzie.
 GEOtop is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __RECOVERY_H__
#define __RECOVERY_H__

#include "rw_maps.h"
#include "struct.geotop.09375.h"
#include "turtle.h"
#include "output.09375.h"
#include "write_ascii.h"

#include <string.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

extern T_INIT *UV;

/* Start variables that have linkage with output.h */
extern double wt0_basin; /*mean intercepted precipitation [mm] in the previous output-basin Dt*/
extern double Ssup;       /*supercial Storage of water in all the basin [mm]*/
extern double Ssub;      /*subsuperficial Storage of water in all the basin [mm]*/
extern double Rout;      /*sum of the output flows from the last output-basin for unit of area[mm]*/
extern double R_G;
extern double S_ch0;     /*wat in the channel at the start of a basin-output step-time*/
extern double S_ch1;     /*wat in the channel at the end of a basin-output step-time*/
extern double Qsub_ch, Qsup_ch ,Q_G; /*averaged output flows*/
extern double SWE_previous;
extern double GWE_previous;
extern double Smelt;   /*Snow melt [mm] during the time interval*/
extern double Ssubl;   /*Snow sublimation [mm] during the time interval*/
extern double Sevap;   /*Snow evaporation [mm] during the time interval*/
extern double Gmelt;   /*Glacier melt [mm] during the time interval*/
extern double Gsubl;   /*Glacier sublimation [mm] during the time interval*/
extern double Gevap;   /*Glacier evaporation [mm] during the time interval*/
extern double Smelt_previous;
extern double Ssubl_previous;
extern double Sevap_previous;
extern double Gmelt_previous;
extern double Gsubl_previous;
extern double Gevap_previous;
extern long isavings;
/* End variables that have linkage with output.h */

unsigned int filecounter;

/** A function to restore a previous state of a simulation. For a given timestamp and on condition that recovery
 *  recovery files for that timestamp exist, the GEOtop structures TIMES, WATER, CHANNEL, SOIL, ENERGY, SNOW and 
 *  GLACIER are overwritten with the values that are stored in the files.
 *
 *  @param timestamp A string representing the timestamp to recover in ISO format (e.g. 20090807T0900)
 *  @param times     A pointer to a TIMES structure, previously allocated
 *  @param wat       A pointer to a WATER structure, previously allocated
 *  @param cnet      A pointer to a CHANNEL structure, previously allocated
 *  @param sl        A pointer to a SOIL structure, previously allocated
 *  @param egy       A pointer to a ENERGY structure, previously allocated
 *  @param snow      A pointer to a SNOW structure, previously allocated
 *  @param glac      A pointer to a GLACIER structure, previously allocated
 */
void recover_simulation(char *timestamp, TIMES *times, WATER *wat, CHANNEL *cnet, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac);


/** A function to write the current state of a simulation. The GEOtop structures TIMES, WATER, CHANNEL, SOIL, 
 *  ENERGY, SNOW and GLACIER are written into files in the recovery folder (./output/rec/).
 *
 *  @param timestamp A string representing the timestamp that will be a suffix for the recovery files 
 *                   in ISO format (e.g. 20090807T0900)
 *  @param times     A pointer to a TIMES structure
 *  @param wat       A pointer to a WATER structure
 *  @param cnet      A pointer to a CHANNEL structure
 *  @param sl        A pointer to a SOIL structure
 *  @param egy       A pointer to a ENERGY structure
 *  @param snow      A pointer to a SNOW structure
 *  @param glac      A pointer to a GLACIER structure
 */
void write_recovery_files(char *timestamp, TIMES *times, WATER *wat, CHANNEL *cnet, SOIL *sl, ENERGY *egy, SNOW *snow, GLACIER *glac);

/** 
 *  @brief A function that restores a DOUBLETENSOR by reading consecutive 2D maps from files
 *  @param _filename  A string representing the beginning of the filename to be recovered
 *  @param timestamp  A string representing a timestamp in ISO format, suffixed to _filename
 *  @param tensor_out A pointer to a pointer of a DOUBLETENSOR
 */
void read_tensor3D(char *_filename, char *timestamp, DOUBLETENSOR **tensor_out);

/** 
 *  @brief A function that restores a SHORT-, LONG- or DOUBLEMATRIX by reading a 2D grass map from a file
 *  @param _filename  A string representing the beginning of the filename to be recovered
 *  @param type       A short representing the type of matrix to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param matrix_out A pointer to a pointer of a matrix; needs to be casted to a (void **);
 */
void read_grass_map(char* _filename, short type, void **matrix_out);

/** 
 *  @brief A function that serves as wrapper for read_grass_map (it builds the filename string)
 *  @param _filename  A string representing the beginning of the filename to be recovered
 *  @param timestamp  A string representing a timestamp in ISO format, suffixed to _filename
 *  @param type       A short representing the type of matrix to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param matrix_out A pointer to a pointer of a matrix; needs to be casted to a (void **);
 */
void read_map2D(char *filename, char *timestamp, short type, void **matrix_out);

/** 
 *  @brief A function that serves as wrapper for write_grassascii (it builds the filename string)
 *  @param _filename  A string representing the beginning of the filename to be recovered
 *  @param timestamp  A string representing a timestamp in ISO format, suffixed to _filename
 *  @param type       A short representing the type of matrix to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param matrix_out A pointer to a pointer of a matrix; needs to be casted to a (void **);
 */
void write_map2D(char *filename, char *timestamp, short type, void *M);

/** 
 *  @brief A function that writes one short value to a file
 *  @param filename  A string representing the complete filename to write to
 *  @param d         A short that shall be written to a file
 */
void write_short_to_file(char* filename, short d);

/** 
 *  @brief A function that writes one long value to a file
 *  @param filename  A string representing the complete filename to write to
 *  @param l         A long int that shall be written to a file
 */
void write_long_to_file(char* filename, long l);

/** 
 *  @brief A function that writes one double value to a file
 *  @param filename  A string representing the complete filename to write to
 *  @param d         A double that shall be written to a file
 */
void write_double_to_file(char* filename, double d);

/** 
 *  @brief A function that reads one short value from a file
 *  @param filename  A string representing the complete filename to read from
 *  @param *s        A pointer to a short that shall be read from the file 
 */
void read_short_from_file(char* filename, short* s);

/** 
 *  @brief A function that reads one long int value from a file
 *  @param filename  A string representing the complete filename to read from
 *  @param *l        A pointer to a long that shall be read from the file 
 */
void read_long_from_file(char* filename, long* l);

/** 
 *  @brief A function that reads one double value from a file
 *  @param filename  A string representing the complete filename to read from
 *  @param *d        A pointer to a double that shall be read from the file 
 */
void read_double_from_file(char* filename, double* d);

/** 
 *  @brief A function that serves as wrapper for read_short_from_file (it builds the filename string)
 *  @param *filename   A string representing the beginning of the filename to be recovered
 *  @param *timestamp  A string representing a timestamp in ISO format, suffixed to filename
 *  @param *s          A pointer to a short that shall be read from the file 
 */
void read_short(char *filename, char *timestamp, short *d);

/** 
 *  @brief A function that serves as wrapper for read_long_from_file (it builds the filename string)
 *  @param *filename   A string representing the beginning of the filename to be recovered
 *  @param *timestamp  A string representing a timestamp in ISO format, suffixed to filename
 *  @param *l          A pointer to a long that shall be read from the file 
 */
void read_long(char *filename, char *timestamp, long *l);

/** 
 *  @brief A function that serves as wrapper for read_double_from_file (it builds the filename string)
 *  @param *filename   A string representing the beginning of the filename to be recovered
 *  @param *timestamp  A string representing a timestamp in ISO format, suffixed to filename
 *  @param *d          A pointer to a double that shall be read from the file 
 */
void read_double(char *filename, char *timestamp, double *d);

/** 
 *  @brief A function that writes a vector to a specified file
 *  @param *filename  A string representing the filename to be recovered
 *  @param type       A short representing the type of vector to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param *vec_in    A void pointer to a VECTOR struct, that will be written to the file
 */
void write_vector_to_file(char* filename, short type, void *vec_in);

/** 
 *  @brief A function that reads a vector (SHORTVECTOR, LONGVECTOR, DOUBLEVECTOR) from a specified file
 *  @param *filename  A string representing the filename to be recovered
 *  @param type       A short representing the type of vector to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param **vec_out  A void pointer to a pointer of a VECTOR struct, that will be recovered; needs to be casted to a (void **);
 */
void read_vector_from_file(char* filename, short type, void **vec_out);

/** 
 *  @brief A function that serves as a wrapper for read_vector_from_file; it constructs the filename with the timestamp
 *  @param *filename  A string representing the filename to be recovered
 *  @param *timestamp A string representing a timestamp in ISO format, suffixed to filename
 *  @param type       A short representing the type of vector to be recovered (0=DOUBLE, 1=LONG, 2=SHORT)
 *  @param **vec_out  A void pointer to a pointer of a VECTOR struct, that will be recovered; needs to be casted to a (void **);
 */
void read_vector1D(char *filename, char *timestamp, short type, void **vector_out);


/**
 * @brief Read all extern declared variables used within output.h from a file called 'output_variables_YYYYMMDDTHHMM'
 * @param timestamp  A string representing a timestamp in ISO format, suffixed to 'output_variables_'
 */
void read_output_variables(char* timestamp);

/**
 * @brief Write all extern declared variables used within output.h into a file called 'output_variables_YYYYMMDDTHHMM'
 * @param timestamp  A string representing a timestamp in ISO format, suffixed to 'output_variables_'
 */
void write_output_variables(char* timestamp);

/**
 * The following helper functions can compare two matrices or tensors and report any differences
 */
void compare_tensor(DOUBLETENSOR *one, DOUBLETENSOR *two);
void compare_matrix(DOUBLEMATRIX *one, DOUBLEMATRIX *two);
void compare_lmatrix(LONGMATRIX *one, LONGMATRIX *two);

/** 
 *  @brief A function that checks whether a file is accessible for reading
 *  @param *filename  A string representing the filename to be checked
 *  @return 0 if file does not exist or is not accessible, 1 otherwise
 */
short file_exists(char *filename);

/** 
 *  @brief A function that serves as a wrapper for file_exists, incrementing 
 *         the variable filecounter if the file exists
 *  @param *filename  A string representing the filename 
 *  @return 0 if file does not exist or is not accessible, 1 otherwise
 */
short check_file_exists(char* filename);


#endif
