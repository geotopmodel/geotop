
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

				
void tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e);

void tridiag2(short a, long r, long c, long nx, double wdi, DOUBLEVECTOR *diag_inf, double wd, DOUBLEVECTOR *diag, double wds, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e);

double norm_2(DOUBLEVECTOR *V, long n);

double norm_inf(DOUBLEVECTOR *V, long n);

double sum(DOUBLEVECTOR *V, long n);

void Cramer_rule(double A, double B, double C, double D, double E, double F, double *x, double *y);

double minimize_merit_function(double res0, double lambda1, double res1, double lambda2, double res2);

double product(DOUBLEVECTOR *a, DOUBLEVECTOR *b);
