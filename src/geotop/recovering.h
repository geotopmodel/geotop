
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

void assign_recovered_map(short old, long n, char *name, Matrix<double> *assign,
                          PAR *par, Matrix<double> *Zdistr);

void assign_recovered_map_vector(short old, long n, char *name, Vector<double> *assign,
                                 Matrix<long> *rc, PAR *par, Matrix<double> *Zdistr);

void assign_recovered_map_long(short old, long n, char *name, Matrix<long> *assign,
                               PAR *par, Matrix<double> *Zdistr);

void assign_recovered_tensor(short old, long n, char *name, DOUBLETENSOR *assign,
                             PAR *par, Matrix<double> *Zdistr);

void assign_recovered_tensor_vector(short old, long n, char *name, Matrix<double> *assign,
                                    Matrix<long> *rc, PAR *par, Matrix<double> *Zdistr);

void assign_recovered_tensor_channel(short old, long n, char *name, Matrix<double> *assign,
                                     Vector<long> *r, Vector<long> *c, Matrix<double> *Zdistr);

void recover_run_averages(short old, Matrix<double> *A, char *name, Matrix<double> *LC,
                          Matrix<long> *rc, PAR *par, long n);

void print_run_averages_for_recover(Matrix<double> *A, char *name, long **j_cont, PAR *par,
                                    long n, long nr, long nc);
