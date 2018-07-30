
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

void i_lrc_cont(Matrix<double> *LC, long ***i, Matrix<long> *lrc, long nl,
                long nr, long nc);

void j_rc_cont(Matrix<double> *LC, long **j, Matrix<long> *rc, long nr, long nc);

void lch3_cont(long **ch3, Matrix<long> *lch, long nl, long nch);

void cont_nonzero_values_matrix2(long *tot, long *totdiag, CHANNEL *cnet,
                                 Matrix<double> *LC, Matrix<long> *lrc, long ***i, long n, long nch, long nl);

void cont_nonzero_values_matrix3(Vector<long> *Lp, Vector<long> *Li,
                                 CHANNEL *cnet, Matrix<double> *LC, Matrix<long> *lrc, long ***i, long n, long nch,
                                 long nl);

