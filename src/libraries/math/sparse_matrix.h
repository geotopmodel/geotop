
/* STATEMENT:

 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 1.225-15 - 20 Jun 2013

 Copyright (c), 2013 - Stefano Endrizzi

 This file is part of Geotop 1.225-15

 Geotop 1.225-15  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

 Geotop 1.225-15  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.

 If you have satisfactorily used the code, please acknowledge the authors.

 */

typedef double (*t_Matrix_element_with_voidp)(long i, DOUBLEVECTOR *eta,
                                              DOUBLEVECTOR *eta0, double dt, void *data);
/****************************************************************************************************/

int get_diagonal(DOUBLEVECTOR *diagonal, DOUBLEVECTOR *x0, double dt,
                 t_Matrix_element_with_voidp Matrix, void *data);
/****************************************************************************************************/

long CG(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x,
        DOUBLEVECTOR *x0, double dt,
        DOUBLEVECTOR *b, t_Matrix_element_with_voidp function, void *data);
/****************************************************************************************************/

int get_upper_diagonal(DOUBLEVECTOR *udiagonal, DOUBLEVECTOR *x0, double dt,
                       t_Matrix_element_with_voidp Matrix, void *data);
/****************************************************************************************************/

long BiCGSTAB(double tol_rel, double tol_min, double tol_max, DOUBLEVECTOR *x,
              DOUBLEVECTOR *x0,
              double dt, DOUBLEVECTOR *b, t_Matrix_element_with_voidp function, void *data);
/****************************************************************************************************/

void product_using_only_lower_diagonal_part(DOUBLEVECTOR *product,
                                            DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);
/****************************************************************************************************/

void product_using_only_upper_diagonal_part(DOUBLEVECTOR *product,
                                            DOUBLEVECTOR *x, LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);
/****************************************************************************************************/

void lower_matrix_product(DOUBLEVECTOR *product, DOUBLEVECTOR *x,
                          LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);
/****************************************************************************************************/

void upper_matrix_product(DOUBLEVECTOR *product, DOUBLEVECTOR *x,
                          LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);

/****************************************************************************************************/

long BiCGSTAB_diag(double tol_rel, double tol_min, double tol_max,
                   DOUBLEVECTOR *x, DOUBLEVECTOR *b,
                   LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
/****************************************************************************************************/

long BiCGSTAB_unpreconditioned(double tol_rel, double tol_min, double tol_max,
                               DOUBLEVECTOR *x, DOUBLEVECTOR *b,
                               LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);
/****************************************************************************************************/

void get_diag_lower_matrix(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag,
                           LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);
/****************************************************************************************************/

void get_diag_upper_matrix(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag,
                           LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);

/****************************************************************************************************/
long BiCGSTAB_lower(double tol_rel, double tol_min, double tol_max,
                    DOUBLEVECTOR *x, DOUBLEVECTOR *b,
                    LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);

/****************************************************************************************************/
long BiCGSTAB_upper(double tol_rel, double tol_min, double tol_max,
                    DOUBLEVECTOR *x, DOUBLEVECTOR *b,
                    LONGVECTOR *Ai, LONGVECTOR *Ap, DOUBLEVECTOR *Ax);

/****************************************************************************************************/

void solve_upper_diagonal_system(DOUBLEVECTOR *x, DOUBLEVECTOR *B,
                                 LONGVECTOR *Ait, LONGVECTOR *Apt, DOUBLEVECTOR *Axt);
/****************************************************************************************************/

void solve_lower_diagonal_system(DOUBLEVECTOR *x, DOUBLEVECTOR *B,
                                 LONGVECTOR *Ait, LONGVECTOR *Apt, DOUBLEVECTOR *Axt);

/****************************************************************************************************/
void solve_SSOR_preconditioning(double omeg, DOUBLEVECTOR *x, DOUBLEVECTOR *B,
                                LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx,
                                LONGVECTOR *Ui, LONGVECTOR *Up, DOUBLEVECTOR *Ux);

/****************************************************************************************************/
long BiCGSTAB_LU_SSOR(double omeg, double tol_rel, double tol_min,
                      double tol_max, DOUBLEVECTOR *x, DOUBLEVECTOR *b,
                      LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx, LONGVECTOR *Ui,
                      LONGVECTOR *Up, DOUBLEVECTOR *Ux);

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/
//FURTHER PERSONALIZATION OF THE FUNCTIONS ABOVE

void product_using_only_strict_lower_diagonal_part(DOUBLEVECTOR *product,
                                                   DOUBLEVECTOR *x, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);

/****************************************************************************************************/
void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(
  DOUBLEVECTOR *product, DOUBLEVECTOR *x, DOUBLEVECTOR *y,
  LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);

/****************************************************************************************************/
void get_diag_strict_lower_matrix_plus_identity_by_vector(DOUBLEVECTOR *diag,
                                                          DOUBLEVECTOR *udiag, DOUBLEVECTOR *y,
                                                          LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);

/****************************************************************************************************/
long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel,
                                                          double tol_min, double tol_max, DOUBLEVECTOR *x,
                                                          DOUBLEVECTOR *b, DOUBLEVECTOR *y, LONGVECTOR *Li, LONGVECTOR *Lp,
                                                          DOUBLEVECTOR *Lx);

/****************************************************************************************************/
void product_matrix_using_lower_part_by_vector_plus_vector(double k,
                                                           DOUBLEVECTOR *out, DOUBLEVECTOR *y, DOUBLEVECTOR *x,
                                                           LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);

/****************************************************************************************************/
void product_using_only_lower_diagonal_part2(DOUBLEVECTOR *product,
                                             DOUBLEVECTOR *x, LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx);