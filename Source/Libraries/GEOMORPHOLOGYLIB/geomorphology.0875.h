/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of part in which you want divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti                                                    */
void sky_view_factor(DOUBLEMATRIX *sky, long N, T_INIT *UV, DOUBLEMATRIX *input, SHORTMATRIX *convess);



void pits_filler_0875(DOUBLEMATRIX *Z0,SHORTMATRIX *land_use);



void nablaquadro_mask(DOUBLEMATRIX *Z0,SHORTMATRIX *curv,DOUBLEVECTOR *U,DOUBLEVECTOR *V);



void nablaquadro(DOUBLEMATRIX *Z0,DOUBLEMATRIX *nabla,DOUBLEVECTOR *U,DOUBLEVECTOR *V);



/* Calculation of the aspect for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - azimuth   matrix with the aspect
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program aspetto of Pegoretti; this subroutine is more indipendent and needs
   less input                                                                   */
void aspect0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *azimuth);



/* Calculation of the mean slopes for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - slopes    matrix with the mean slopes
   Subroutine created by Davide Tamanini (June 2003)                                */
void slopes0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *slopes);



/*Computation of the thickness of hydrological active soil (bedrock depth) using
  the method of Dietrich:
   Input:  - soil_parameters  par of soil (in particular that of Dietrich)
           - UV               struct with vector of pixel dimensions and novalue
           - Z0               matrix with elevation (DTM)
   Output: - h                matrix with the depth of bedrock [mm]
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   subroutine profondita of Pegoretti; this subroutine is more indipendent
   and needs less input; it was modified by Davide Tamanini (April 2004)
   to eliminate the matrix of dranaige direction as input.                      */
void soil_depth(double min_h, double h_crit, double P, double M, double prs, double k , T_INIT *UV, DOUBLEMATRIX *Z0, DOUBLEMATRIX *h);



/* Calculation of area considering the slopes for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - area      matrix with the mean slopes
   Subroutine created by Davide Tamanini (June 2003)                                */
void area0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *area);



/* Computation of the Drainage direction following this convention:
   0 for no draining pixel; from 1 to 8 starting from east in toward counterclockwise;
   9 for novalue pixel; 10 for outlet pixel; 11 for lake and sea; 12 pit pixel
   Inputs   - Z0               matrix with elevation (DTM)
            - land_use         matrix with land use
            - UV               struct with vector of pixel dimensions and novalue
   Outputs:	- directions       matrix with Drainage Direction                         */
void DrainageDirections0875(DOUBLEMATRIX *elevations,SHORTMATRIX *land_use,T_INIT *UV,SHORTMATRIX *directions);



//presa uguale da geomorphology099
void gradients(DOUBLEMATRIX *Z0,SHORTMATRIX *directions,DOUBLEMATRIX *dr_gradient,T_INIT *UV);



/* Calculation of channel network: this subroutine is a modification of select_hillslopes
   created by Davide Tamanini (August 2003)                                                */
void select_hillslopes_mod(LONGMATRIX *ca,DOUBLEMATRIX *dr_gradient,DOUBLEMATRIX *curv,SHORTMATRIX *m,double threshold,DOUBLEVECTOR *pixelsize);



/*Computation of the trasversal slope of the channel pixel (created by Davide 30-8-2003):*/
void channel_lateral_slope(DOUBLEMATRIX *Z0,SHORTMATRIX *DD,T_INIT *UV,DOUBLEMATRIX *i_ch);



void distance_from_channel(DOUBLEMATRIX *dist, SHORTMATRIX *DD, SHORTMATRIX *ST);

void distance_from_channel2(DOUBLEMATRIX *dist, SHORTMATRIX *ST, LONGVECTOR *rch, LONGVECTOR *cch);

void set_boundary_condition(DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, short code, SHORTMATRIX *pixel_type, double novalue);

