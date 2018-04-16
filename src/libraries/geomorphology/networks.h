#define MAXELEVATION 10000


/**

Name:  tca

Synopsis:  void tca(SHORTMATRIX *,LONGMATRIX *);


Version:  0.9

Description: It calculates the total contributing area for a given point
in a network

Authors & date: Riccardo Rigon, November 1997


Inputs: 1)The matrix of flowing directions; 2) The matrix that will contain
te values of contributing area

Returns: void

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


Notes: The function can be made much more faster using an algorithm
that uses the branching structure of the networks.


*/

void tca(SHORTMATRIX *,LONGMATRIX *);

/**

Name: outletdistance

Synopsis: void outletdistance(SHORTMATRIX *m,DOUBLEMATRIX *dist,DOUBLEVECTOR *U);


Version: 0.9

Description: It returns the distance from the outlet or outlets of the given DEM. Diagonal  flowing directions are given distance according
to the Pitagora theorem.  Outlet pixel is given null distance. Units are pixels.

Authors & Date: Riccardo Rigon, January 1998


Inputs: 1) the flowing direction matrix; 2) the matrix that will contain the distances; 3) The vector of novalues for the flowing directions

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997

*/

void outletdistance(SHORTMATRIX *m,DOUBLEMATRIX *dist,DOUBLEVECTOR *U);


/**

Name: topological_outletdistance

Synopsis: void topological_outletdistance(SHORTMATRIX *m,LONGMATRIX *dist)

Version: 0.9

Inputs: 1) the flowing direction matrix; 2) the matrix that will contain the distances; 3) The vector of novalues for the flowing directions

Description: It returns the distance from the outlet or outlets of the given DEM. Units are pixels.
Diagonal  flowing directions are given distance 1.  Outlet pixel is given null distance.


Authors & Date: Riccardo Rigon, January 1998

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

void topological_outletdistance(SHORTMATRIX *m,LONGMATRIX *dist);

/**

Name:  tplgicl_hc_outletdistance

Synopsis: void topological_hillslopes_channels_outletdistance(SHORTMATRIX *m,SHORTMATRIX* hillslopes,
               LONGMATRIX *dist)

Version: 0.1

Description:  It returns the topological distance from the outlet .
Hillslope pixels are marked with a  negative number whose absolute value is the distance
from channels mesured along the pathway identified by steepest descent.    Outlet pixel has null distance.

Authors & Date: Riccardo Rigon, January 1998

Inputs: 1) the flowing direction matrix; 2) the matrix that distinguish hillslopes from channels ;3) the matrix that will contain the distances

See Also: outlet_distance, topological_outletdistance

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

void topological_hillslopes_channels_outletdistance(SHORTMATRIX *m,
                                                    SHORTMATRIX *hillslopes,
                                                    LONGMATRIX *dist);

/**

Name: h_c_outletdistance

Synopsis: void hillslopes_channels_outletdistance(SHORTMATRIX *m,SHORTMATRIX* hillslopes,
              DOUBLEMATRIX *dist,DOUBLEVECTOR *U)


Version: 0.9

Description:   It returns the distance from the outlet .
Hillslope pixels are marked with a  negative number whose absolute value is the distance
from channels mesured along the pathway identified by steepest descent.   Outlet pixel has null distance.

Authors & Date: Riccardo Rigon, January 1998

Inputs: 1) the flowing direction matrix; 2) the matrix that marks hillslopes; 3) the matrix that will contain the distances; 4) The vector of novalues for the flowing directions

See Also: outlet_distance, topological_outletdistance

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

void hillslopes_channels_outletdistance(SHORTMATRIX *m,
                                        SHORTMATRIX *hillslopes, DOUBLEMATRIX *dist,DOUBLEVECTOR *U);
/* old veriosn with U float */
void hillslopes_channels_outletdistance_f(SHORTMATRIX *m,
                                          SHORTMATRIX *hillslopes, DOUBLEMATRIX *dist,FLOATVECTOR *U);
/**

Name: select_channel

Synopsis: void select_channels(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,DOUBLEMATRIX *curv,
                       SHORTMATRIX *m,double threshold);


Version: 0.9

Description:   It  uses curvatures and contributing areas to estract channels from drainage directions. i.e. those points that
have contributing area larger than threshold and have positive curvature are marked as channels.

Authors & Date: Riccardo Rigon, January 1998

See Also:  select_hillslope

Inputs: 1) the contributing area matrix; 2) the drainage directions; 3) the curvature matrix;
4) the matrix that will contain the channels pixels marked; 5) a threshold value


FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

void select_channels(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,
                     DOUBLEMATRIX *curv,
                     SHORTMATRIX *m,double threshold);


/**

Name: select_hillslopes

Synopsis:  void select_hillslopes(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,DOUBLEMATRIX *curv,
                       SHORTMATRIX *m, double threshold);


Version: 0.9

Description:   It  uses curvatures and contributing areas to estract channels from drainage directions.
  i.e. those points that  have contributing area larger than threshold and have negative curvature
  are marked as hillslope.

Authors & Date: Riccardo Rigon, January 1998

See Also:  select_channels

Inputs: 1) the contributing area matrix; 2) the drainage directions; 3) the curvature matrix;
4) the matrix that will contain the channels pixels marked; 5) a threshold value


FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,

References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997



*/

void select_hillslopes(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,
                       DOUBLEMATRIX *curv,
                       SHORTMATRIX *m, double threshold);
/* old version with float for GEOTOP 075 */
void select_hillslopes_f(LONGMATRIX *ca,FLOATMATRIX *drainagedirections,
                         FLOATMATRIX *curv,
                         SHORTMATRIX *m, double threshold);

/**

Name:  randomflow

Synopsis:  void randomflow(long ,long , SHORTMATRIX *);


Version:  0.9

Description: given the matrix of flowing directions and a position
in the matrix it returns a new flowing direction selected at random

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) 2) the position; 3) The pointer to the flowing direction matrix

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997




*/
/* commented because the Warning: variable 'work' is not initialized before being used
void randomflow(long ,long , SHORTMATRIX *);

/**

Name:  crossQ

Synopsis:  long crossQ(long  ,long ,long ,long ,SHORTMATRIX *);



Version:  0.9

Description: any new flowing direction must be consistent in the sense
that it cannot across older flower directions. crossQ checks if the above case
happens.

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) 2) the position; 3) 4) the coordinates chosen as new flowing
direction. 5)  The pointer to the flowing direction matrix

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

long crossQ(long,long,long,long,SHORTMATRIX *);


/**

Name:  coords2int

Synopsis:  short coords2int(long ,long ,long ,long );


Version:  0.9

Description: A flowing direction is saved as a number from 0 to 9.
coords2 int performs the transformation among the coordinate of the
point, the flowing directions and the above number.

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) 2) the position; 3) 4) the flowing position;

Returns: a number between 0 and 9 indicating the flowing direction



Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997

*/

short coords2int(long,long,long,long );

/**

Name:  go_dowstream

Synopsis:  void go_downstream(long *,short );


Version:  0.9

Description: replace the coordinates of a point in a river basin
with the coordinates of the point downstream

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) the position; 2) the flowing direction (coded); 3) the number
of columns in matrix (It is not used if the constant CILINDRICAL is not set to 1)

Returns: void



Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997



*/

void go_downstream(long *,short,long );
void go_downstream_2(long *,short,long );





/**

Name:  go_upstream_a

Synopsis: long go_upstream_a(double th,long *p, long *kk,SHORTMATRIX *m,LONGMATRIX *tca,DOUBLEMATRIX *l);




Version:

Description: go_upstream_a definisce il numero di pixel che drenano nel pixel
       definito da *p. La routine definisce con p di uscita il pixel di drenaggio
       che rappresenta il ramo princiapale. Nel caso vi siano piu' pixel drenanti
       quello che definisce il ramo princiaple viene definito usando l'area cumulata e la lunghezza
       dando prevalenza alla prima


Authors & date: Riccardo Rigon, Marco Pegoretti Luglio 1999

Inputs: 1)matrix of the directions
           2)matrix of the aree
           3)matrix of the l
           4)vettore che indica la posizione del pixel
           5)kk indica la direzione in cui si trova il pixel di drenaggio principale
           6)threshold of the aree


Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

long go_upstream_a(double th,long *p, long *kk,SHORTMATRIX *m,LONGMATRIX *tca,
                   DOUBLEMATRIX *l);

short go_upstream_b(long *p,SHORTMATRIX *,LONGMATRIX *,DOUBLEMATRIX *);

void initialize_longmatrix(LONGMATRIX *,long );

void initialize_flowdirections(SHORTMATRIX * );

LONGPAIR *initialize_flowdirections_with_outlets(SHORTMATRIX *m);

/**

Name:  energyexpenditure, weighted_energyexpenditure,r_energyexpenditure,
r_weighted_energyexpenditure

Synopsis:

double energyexpenditure(LONGMATRIX *,double ,long );
double weighted_energyexpenditure(SHORTMATRIX *,LONGMATRIX *,double ,long );
double r_energyexpenditure(DOUBLEMATRIX *tca,double ex,long th)



Version:  0.9

Description: returns the energy expenditure in a river basin in pixel units

Authors & date: Riccardo Rigon, November 1997


Inputs: 1) contributing area matrix;

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/
double energyexpenditure(LONGMATRIX *,double,long );

double weighted_energyexpenditure(SHORTMATRIX *,LONGMATRIX *,double,long );

double r_energyexpenditure(DOUBLEMATRIX *tca,double ex,double th);

double r_weighted_energyexpenditure(SHORTMATRIX *,DOUBLEMATRIX *,double,
                                    double );

/**

Name:  topology

Synopsis:  LONGPAIR * topology(void );



Version:  0.9

Description: Used by metropolis to know which topology is being to be used in simulation. So far only the D8 topology is implemented
but it should be simple to generalize the choice to D4.

Authors & date: Riccardo Rigon, November 1997

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/
LONGPAIR *topology(void );

/**

Name: tcamax

Synopsis: long tcamax(SHORTMATRIX *,LONGMATRIX *,DOUBLEMATRIX *,long *,long ,double );

Version: 1.0

Description: used by hacklengths.
tcamax scans the neighborhood of a  given point, P, and returns 1 (true) if one
of the points that drain into P has greater tca than the  specified. In  the case
of two points with the same tca, the more distant from the divides is chosen.
It returns 0 otherwise.




Inputs:
1) The pointer to the matrix of flowing directions;
2) The pointer to matrix of weigths;
3) The pointer to the matrix of distances;
4) The pointer to a vector (in normal C sense) containing the point P;
5) The value of the maximum of contributing area;
6) The distance of the point P from the divides


Authors & date: Riccardo Rigon, November 1997

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


See Also: hacklenghs, sourcesq


*/

long tcamax(SHORTMATRIX *,LONGMATRIX *,DOUBLEMATRIX *,long *,long,double );


/**



Name: sourcesq

Synopsis: long sourcesq(SHORTMATRIX *,long *);

Version: 1.0

Description: It return 1 if the point is a source, 0 otherwise


Authors & date: Riccardo Rigon, November 1997


Inputs:
1) The pointer to the matrix of flowing directions;
2) The pointer to a vector (in normal C sense) containing the point P;


Authors & date: Riccardo Rigon, November 1997

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

long sourcesq(SHORTMATRIX *,long *);

/**



Name: linkq, conditioned_linkq

Synopsis:

long linkq(SHORTMATRIX *m,long *flow);
long conditioned_linkq(SHORTMATRIX *m,DOUBLEMATRIX *ca,long *flow,double th)


Version: 1.0

Description: It return 1 if the point is a link, 0 otherwise. conditioted_linkq
tests also if the ca field value in the chosen point is greater than a threshold
value


Authors & date: Riccardo Rigon, October 1998


Inputs:
For linkq:
1) The pointer to the matrix of flowing directions;
2) The pointer to a vector (in normal C sense) containing the point P;
For conditioned_linkq:
1) The pointer to the matrix of flowing directions;
2) The pointer to the matrix of the testing fields;
3) The pointer to a vector (in normal C sense) containing the point P;
4) A threshold value

Authors & date: Riccardo Rigon, November 1997

Examples:  metropolis.c,  eden.c

FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997

See Also:   sourcesq




*/
long linkq(SHORTMATRIX *m,long *flow);
long conditioned_linkq(SHORTMATRIX *m,DOUBLEMATRIX *ca,long *flow,double th);

/**



Name: downstream

Synopsis: void downstream(SHORTMATRIX *flow,DOUBLEMATRIX *quantity,DOUBLEMATRIX *laggedquantity,long lag,double dx, double dy);

Version: 0.8

Description: downstream fills 'laggedquantity' with   the values of 'quantity' contained 'lag' times downstream
for each point in the original quantity. The result can be used to calculate covariances or variograms of the same
quantity or related quantities with lags measured in the steepest descent direction. Diagonl directions are weighted
according to the square root of 2.

Authors & date: Riccardo Rigon, February 1998


Inputs:
1) the pointer to the matrix of flowing directions;
2) the pointer to one the matrix of which one wants to calculate the autocorrelation;
3) the pointer to the matrix containing the lagged values;
4) the lag;
5) the length of a pixel in the x direction;
6) the length of a pixel in the y direction;




FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997




*/

void downstream(SHORTMATRIX *,DOUBLEMATRIX *,DOUBLEMATRIX *,long lag,double,
                double,double);

double network_dowstreamcorrelation(DOUBLEMATRIX *,DOUBLEMATRIX *,double,
                                    double,double,double );


/**



Name: drainagedirections



Version: 1.0

Synopsis: void drainagedirections(DOUBLEMATRIX *,SHORTMATRIX *, DOUBLEVECTOR* U,DOUBLEVECTOR* V);


Description: DrainageDirections find the steepest descent mark it as drainage direction


Authors & date: Riccardo Rigon, Paolo Verardo, October 1997



Inputs: 1) The matrix of elevations; 2) The empty matrix of directions; 3)
a vector of double containing the pixel sizes.


FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

void drainagedirections(DOUBLEMATRIX *elevations,SHORTMATRIX *directions,
                        DOUBLEVECTOR *U,DOUBLEVECTOR *V);




/**



Name: drainagedirections_modify



Version: 1.0

Synopsis: void drainagedirections_modify(DOUBLEMATRIX *,SHORTMATRIX *, DOUBLEVECTOR* U,DOUBLEVECTOR* V
                                                            SHORTMATRIX *);


Description: DrainageDirections find the steepest descent and marks it as drainage direction


Authors & date: Riccardo Rigon, Paolo Verardo, October 1997



Inputs: 1) The matrix of elevations; 2) The empty matrix of directions; 3)
a vector of double containing the pixel sizes.





FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997

*/

void drainagedirections_modify(DOUBLEMATRIX *,SHORTMATRIX *, DOUBLEVECTOR *U,
                               DOUBLEVECTOR *V,
                               SHORTMATRIX *);

/* old version with U, V floating point */
void drainagedirections_modify_f(DOUBLEMATRIX *,SHORTMATRIX *, FLOATVECTOR *U,
                                 FLOATVECTOR *V,
                                 SHORTMATRIX *);

/**



Name: is_ontheborder



Version: 1.0

Synopsis: short is_ontheborder(DOUBLEMATRIX *elevations,DOUBLEVECTOR *V,long i, long j);


Description: is_ontheborder return 1 if the point is a basin border point
0 otherwise

Authors & date: Riccardo Rigon, May 1999



Inputs: 1) The matrix of elevations; 2) the novalues vector; 3)4)
the row and column of the point





Returns: 1 is succesfull 0 otherwise



FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References: I. Rodriguez_Iturbe and A. Rinaldo, Fractal River Basins, CUP 1997


*/

short is_ontheborder(DOUBLEMATRIX *elevations,DOUBLEVECTOR *V,long i, long j);
/* old version with U, V floating point */
short is_ontheborder_f(DOUBLEMATRIX *elevations,FLOATVECTOR *V,long i,
                       long j);
/**



Name: sum_downstream



Version: 1.0

Synopsis: void sum_downstream(SHORTMATRIX *flow,DOUBLEMATRIX *la,DOUBLEMATRIX *dist);


Description: sum_downstream takes a quantity and accumulate it downstream fllowing the steepest descent.

Authors & date: Riccardo Rigon, May 1999



Inputs: 1) The matrix of flow; 2) the matrix containing the value to be summed; 3)4
the output matrix



FILE:  LIBRARIES/GEOMORPHOLOGYLIB/networks.c,   LIBRARIES/GEOMORPHOLOGYLIB/networks.h,


References:

Banavar et al.,  Allometry of River Networks,  submitted for pubblication, 2000


*/


void sum_downstream(SHORTMATRIX *flow,DOUBLEMATRIX *la,DOUBLEMATRIX *dist);

void dem_array_check(DOUBLEVECTOR *T, DOUBLEVECTOR *U );

