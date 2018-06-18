/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti                                                    */
void sky_view_factor(DOUBLEMATRIX *sky, long N, T_INIT *UV,
                     DOUBLEMATRIX *input, SHORTMATRIX *convess, long novalue);


void nablaquadro_mask(DOUBLEMATRIX *Z0,SHORTMATRIX *curv,Vector<double>* U,
                      Vector<double>* V);