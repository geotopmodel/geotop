#define RATIO 1000

/**
 
Name: vectorize_doublematrix

Version: 1.0

Synopsis: DOUBLEVECTOR *vectorize_doublematrix(DOUBLEMATRIX *input)


Description: takes a double matrix and transforms it to a vector. Actually
this is easy done cause the data of the matrix  are stored in adjacent addresses.
Thus, the function allocat the vector structure; put the address of the vector 
elements to the first of the matrix elements. Free the memory used by the matrix structure
(but not that used by data). 

Inputs:  the input matrix

Return: the pointer to the doublevector containing now the data.

Examples:  APPLICATIONS/DATA_MANIPULATION/coupledfield_moments.c

Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

*/

DOUBLEVECTOR *vectorize_doublematrix(DOUBLEMATRIX *input);
FLOATVECTOR *vectorize_floatmatrix(FLOATMATRIX *input);
LONGVECTOR *vectorize_longmatrix(LONGMATRIX *input);
SHORTVECTOR *vectorize_shortmatrix(SHORTMATRIX *input);

/**
 
Name: sortreal

Version: 1.0

Synopsis: void sortreal(DOUBLEVECTOR *ra)


Description: sort a vector of double. The routine is very similar to the heap sort
presented in Numerical Recipes to which one could refer for reference.

Inputs:  the input vector


Examples:  APPLICATIONS/DATA_MANIPULATION/coupledfield_moments.c

Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c


*/

void sortreal(DOUBLEVECTOR *);

/**
 
Name: sort2realvectors
 
Version: 1.0

Synopsis: void sort2realvectors(DOUBLEVECTOR *,DOUBLEVECTOR *);

Description: sort a vector of double and accorsingly a second
vector. The routine is very similar to the heap sort
presented in Numerical Recipes to which one could refer for reference.

Inputs:  the input vector

Examples:  APPLICATIONS/DATA_MANIPULATION/coupledfield_moments.c


Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c
 

*/

void sort2realvectors(DOUBLEVECTOR *,DOUBLEVECTOR *);
void  sort2floatvectors(FLOATVECTOR *ra,FLOATVECTOR *rb);
void  sort2vectors(LONGVECTOR *,FLOATVECTOR *);
/**
 
Name: realpair_intodoublematrix, xyz_into_doublematrix

Synopsis: void realpair_into_doublematrix(REALPAIR * head,DOUBLEMATRIX *indx )
      void xyz_into_doublematrix(XYZ * head,DOUBLEMATRIX *indx )


Description: takes a linked list (of two or three elements) and transform
it into a matrix of double of the appropriate dimensions

Inputs: 1)  the linked list head; 2) the doublematrix pointer



Examples: hystogram.c 

See Also: 

Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c
 

*/




void realpair_into_doublematrix(REALPAIR *,DOUBLEMATRIX *);
void xyz_into_doublematrix(XYZ *,DOUBLEMATRIX *);

/**
 
Name:  simplehystogram, loghystogram

Version: 1.0

Synopsis:
DOUBLEMATRIX *simplehystogram(DOUBLEVECTOR *U,long N,long mn), 
DOUBLEMATRIX *loghystogram(DOUBLEVECTOR *U,long N,long mn, long base)


Description: Doing the hystogram of a set of data is more art than science. 
It requires to discretize the domain of the sample data in  intervals suitable 
for size and shape. Each interval is usually called 'bin'. 
It seems  obvious that the bin have to be of equal size. But it is not. 
Suppose for instance that the data collects the number of points sharing
the same number of upstream sources. Exploring the river in the downstream direction one can see
that this number makes jumps according to a power law. Also, in a river of infinite size,
the values of the possible number of sources span the entire set of integers
numbers. Nevertheless a real river behave like a finite sample of the infinite
one and the domain of the sample will cover only a spotted subset of the integers
where spots (i.e. empty intervals) become larger and larger with size.
This means that keeping fixed the size of the bin most of the bins  will be empty.
A recipe is to take then ( an exponentially) variable size bins. 

Inputs: 1) a sorted vector containing the data to be binned; 2) the number of bins
(instead of fixing the size of each bin is usually convenient to select a fixed number 
of bins ). If the number of bins is set to 0, is simply counted the number of elements with the
same abscissa; 3) the minimum number of elements required in each bin (this implies
that if a bin has not enough elements the numbers of elements of two adjacent bins are summed).
In the case of the exponential binning, the j-th bin extends from base^(j*delta) 
to base^((j+1)*delta) where base is the fourth input field and delta the extension of 
each bin in the logarithm  axis.


Return: a matrix of double: the first column contains the number of elements
in each mean, the second coulumns the mean abscissa of the data in the bin, the third
the highest limit of each bin interval 

Examples: hystogram.c 


Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Notes: In the case of an exponential binning, the mean abscissa
of the data is calculated with the aritmetic mean (it could be calculated as
the geometric mean instead. Because  novalues is not specified novalues are not dropped
from the data vector.


*/


DOUBLEMATRIX *simplehystogram(DOUBLEVECTOR *,long ,long);
DOUBLEMATRIX *exponentialhystogram(DOUBLEVECTOR *,long ,long , long);

/**
 
Name: initialize_

Version: 1.0

Synopsis: 
void initialize_doublevector(DOUBLEVECTOR *,double );
void initialize_longmatrix(LONGMATRIX *, long);


Description:  It initialize the matrix or vector with  the specified value

Inputs:  1) The pointer to the structure to initalize; 2) the value used for initialization

Examples:  Variogram.c


Authors & Date: Riccardo Rigon, October 1997. 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c


*/


void initialize_longvector(LONGVECTOR *, long);
void initialize_shortvector(SHORTVECTOR *,short );
void initialize_intvector(INTVECTOR *,int );
void initialize_floatvector(FLOATVECTOR *,float );
void initialize_doublevector(DOUBLEVECTOR *,double );

void initialize_longmatrix(LONGMATRIX *, long);
void initialize_shortmatrix(SHORTMATRIX *,short );
void initialize_intmatrix(INTMATRIX *,int );
void initialize_floatmatrix(FLOATMATRIX *,float );
void initialize_doublematrix(DOUBLEMATRIX *,double );


/**



Name: copy__matrix

Synopsis: 

void copy_longmatrix(LONGMATRIX *,LONGMATRIX *);
void copy_shortmatrix(SHORTMATRIX *,SHORTMATRIX *);
void copy_intmatrix(INTMATRIX *,INTMATRIX *);
void copy_floatmatrix(FLOATMATRIX *,FLOATMATRIX *);
void copy_doublematrix(DOUBLEMATRIX *,DOUBLEMATRIX *);
void copy_shortvector(SHORTVECTOR *,SHORTVECTOR *);
void copy_intvector(INTVECTOR *,INTVECTOR *);
void copy_longvector(LONGVECTOR *,LONGVECTOR *);
void copy_floatvector(FLOATVECTOR *,FLOATVECTOR *);
void copy_doublevector(DOUBLEVECTOR *,DOUBLEVECTOR *);

Version: 0.8

Description: copy a matrix into a  new one

Authors & date: Riccardo Rigon, January 1998

Inputs: 1- the destination matrix; 2- the origin matrix


FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c




*/




void copy_longmatrix(LONGMATRIX *,LONGMATRIX *);
void copy_shortmatrix(SHORTMATRIX *,SHORTMATRIX *);
void copy_intmatrix(INTMATRIX *,INTMATRIX *);
void copy_floatmatrix(FLOATMATRIX *,FLOATMATRIX *);
void copy_doublematrix(DOUBLEMATRIX *,DOUBLEMATRIX *);
void copy_shortvector(SHORTVECTOR *,SHORTVECTOR *);
void copy_intvector(INTVECTOR *,INTVECTOR *);
void copy_longvector(LONGVECTOR *,LONGVECTOR *);
void copy_floatvector(FLOATVECTOR *,FLOATVECTOR *);
void copy_doublevector(DOUBLEVECTOR *,DOUBLEVECTOR *);
void add_doublevector(DOUBLEVECTOR *small, DOUBLEVECTOR *big);

/**



Name:  _element_multiplication

Synopsis: 
void shortvector_element_multiplication(SHORTVECTOR* ,SHORTVECTOR *);
void floatvector_element_multiplication(FLOATVECTOR* ,FLOATVECTOR *);
void doublevector_element_multiplication(DOUBLEVECTOR* ,DOUBLEVECTOR *);
void longvector_element_multiplication(LONGVECTOR* ,LONGVECTOR *);

void shortmatrix_element_multiplication(SHORTMATRIX* ,SHORTMATRIX *);
void floatmatrix_element_multiplication(FLOATMATRIX* ,FLOATMATRIX *);
void doublematrix_element_multiplication(DOUBLEMATRIX* ,DOUBLEMATRIX *);
void longmatrix_element_multiplication(LONGMATRIX* ,LONGMATRIX *);


Version: 0.1

Description: multiplies each element of the first matrix (vector) 
with the correspondent element of the second. The result overwrite the first matrix

Authors & History: Riccardo Rigon, February 1998


Inputs: the two matrixes or vectors to multiply



References: 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Notes: No control is made on the maximum or minimum value
and calculation can  go out of range

*/
void shortvector_element_multiplication(SHORTVECTOR* ,SHORTVECTOR *);
void floatvector_element_multiplication(FLOATVECTOR* ,FLOATVECTOR *);
void doublevector_element_multiplication(DOUBLEVECTOR* ,DOUBLEVECTOR *);
void longvector_element_multiplication(LONGVECTOR* ,LONGVECTOR *);

void shortmatrix_element_multiplication(SHORTMATRIX* ,SHORTMATRIX *);
void floatmatrix_element_multiplication(FLOATMATRIX* ,FLOATMATRIX *);
void doublematrix_element_multiplication(DOUBLEMATRIX* ,DOUBLEMATRIX *);
void longmatrix_element_multiplication(LONGMATRIX* ,LONGMATRIX *);

/**
 
Name: split, exponentialsplit
 
Synopsis: 

DOUBLEBIN *split(DOUBLEVECTOR *tobesplitted,long N ,double novalue);
DOUBLEBIN *esponentialsplit(DOUBLEVECTOR *tobesplitted,long N ,double base,double novalue);

Version: 0.9

Description: it takes a vector of double and split it  in N parts. 
If N <1 each part just contains the elements of the vectors that have the same values, if N >2
first, the range of data is subdivided in N-1 intervals of equal size from the minimum value (say: min)  to the maximunm value 
(say: max) in the vector escluding those elements that are marked as 'novalue' - a novalue must be either greater or less 
than the effective elements values - . Secondly lists the elements that lie in a bin from min-0.5*delta to max+0.5*delta, for
a total number of N bins. The program works  setting a proper set of pointer to transform the vector in a bin of double (see turtle.h) 
the pointer to the vector is eventually deallocated and substituted (you do not have to deallocate it anymore ) 
with the pointer to  a bin. 
exponentialsplit works exactly the same as split but the bins are of equal size in the logarithmic space, meaning that delta is
exponentially varying  with increasing values. One can chose the logarithm base in which to work.

Inputs: 
split: 1) The vector to be splitted; 2) The number of bins; 3) the novalue
exponentialsplit: 1) The vector to be splitted; 2) The number of bins; 4) The base of the logarithm; 5) the novale

Authors & Hystory: Riccardo Rigon, Paolo D'Odorico, November 1997, February 1998 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Return: a pointer to the DOUBLEBIN

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Needs: LIBRARIES/BASICMATHSTAT/statistics.c, LIBRARIES/BASICMATHSTAT/t_statistics.h


References: Fractal River Basins, I. Rodriguez-Iturbe and A. Rinaldo, C.U.P., 1997

Notes: Probably there is not a correct cleaning on memory ... when trasforming
the vector into the bin the novalue part of the vector remain kind of hidden

*/

DOUBLEBIN *split(DOUBLEVECTOR *,long ,FLOATVECTOR  *);
DOUBLEBIN *esponentialsplit(DOUBLEVECTOR *,long ,double ,FLOATVECTOR  *);


/**
 
Name: split2realvectors, exponentialspliterealvectors

Version: 0.9
 
Synopsis: 

void	split2realvectors(DOUBLEVECTOR *FIRST,DOUBLEVECTOR *SECOND,DOUBLEBIN *ONE,DOUBLEBIN *TWO,long N,double novalue);
void	esponentialsplit2realvectors(DOUBLEVECTOR *FIRST,DOUBLEVECTOR *SECOND,DOUBLEBIN *ONE,DOUBLEBIN *TWO,long N,double base, double  novalue);



Description: They works as split or exponentialsplit on the first of two real vectors. The second one is
splitted in the same position the first one is splitted. 

Inputs: 
split2realvectors: 1) The first vector to be splitted; 2) The second vector to be splitted; 
3) The pointer to the first resulting bin; 4) The pointer to the second resulting bin; 5)  The number of bins; 6) the novalue.
exponentialsplit2realvectors: 1) The first vector to be splitted; 2) The second vector to be splitted; 3) The pointer to the first resulting 
bin; 4) The pointer to the second resulting bin; 5)  The number of bins;6) the base of the logarithm; 7) the novalue.

Authors & Hystory: Riccardo Rigon, Paolo D'Odorico, November 1997, February 1998 

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Needs: LIBRARIES/BASICMATHSTAT/statistics.c, LIBRARIES/BASICMATHSTAT/t_statistics.h

Examples: coupledfield_moments.c 

See Also: hystogram, coupled_field_moments.c

References: Fractal River Basins, I. Rodriguez-Iturbe and A. Rinaldo, C.U.P., 1997

NOTES: Probably there is not a correct cleaning of memory ... when trasforming
the vector into the bin the novalue part of the vector remains hidden. Instead of using a longpair
list, it can be used a simpler type of linked list with only one value stored.
Low testing.

*/

double	split2realvectors(DOUBLEVECTOR *,DOUBLEVECTOR *,DOUBLEBIN *,DOUBLEBIN *,long,long,FLOATVECTOR *);
double    esponentialsplit2realvectors(DOUBLEVECTOR *,DOUBLEVECTOR *, DOUBLEBIN* ,DOUBLEBIN* ,long ,long, double ,FLOATVECTOR * );

/**-----------------------------------------------------------------------*/

/**
 
Name: clean_floatmatrix
 
Synopsis: void clean_floatmatrix(FLOATMATRIX *iv,FLOATMATRIX *ov,FLOATVECTOR *U,FLOATVECTOR *V);


Description: sets to novalue1 those points in the first matrix that are marked as novalue2 in
the second. A reasonable novalue is either smaller or larger than any other value in the matrix.
In the first case, the first element of a novalue vector is <0. In the secon case it is >0. It is 0 if
for some reason the novalue is some value value in between the data values. 

Inputs: 1) the pointer to the first matrix; 2) the pointer to the second matrix; 3)
the novalue of the first matrix; 4) the novalue of the second matrix;

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c


Authors & Date: Riccardo Rigon, March 1, 1998


*/

void clean_floatmatrix(FLOATMATRIX *iv,FLOATMATRIX *ov,FLOATVECTOR *U,FLOATVECTOR *V);

/**
 
Name: shrink_doublematrix
 
Synopsis: DOUBLEMATRIX *shrink_doublematrix(DOUBLEMATRIX *data,FLOATVECTOR *Q);


Description: eliminates from the first matrix (data) those elements set to novalues
and return a smaller matrix

Inputs: 1) the pointer to the data matrix;  2)
the novalue of the first matrix; 

Returns: the new "shrinked" matrix

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Examples: simple_interpolation.c

Authors & Date: Riccardo Rigon, March 1, 1998

References:

Notes: Has been thought for matrixes with 2 columns

*/

DOUBLEMATRIX *shrink_doublematrix(DOUBLEMATRIX *data,FLOATVECTOR *Q);

/**
 
Name: interpolating_function
 
Synopsis: DOUBLEMATRIX *interpolating_function(DOUBLEMATRIX *cleandata);


Description: takes a set of X Y odered pairs and  produces the coefficient of the
line passing trought two  sdiacent pairs. These are stored in a matrix of double

Inputs: 1) the pointer to the data matrix to be interpolated; 

Return: a pointer to the matrix of the angular coefficient and intercept

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Examples: simple_interpolation.c

See Also: interpolate

Authors & Date: Riccardo Rigon, February, 1999

References:  S. Wolfram, Mathematica 3.0, Cambridge University Press


*/

DOUBLEMATRIX *interpolating_function(DOUBLEMATRIX *cleandata);

/**
Name: interpolate
 
Synopsis: double interpolate(double x,DOUBLEMATRIX *cleandata,DOUBLEMATRIX* W);

Description: linearly interpolate the values of a function at x based onthe data
in cleandata and the regression coefficient in w

Inputs: 1) the point where to calculate the function; 2) the data set that serve as basis
for the interpolation; 3) a pointer to the matrix containing the regressions par 

Return: the interpolated value

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Examples: simple_interpolation.c

See Also: interpolating_function

Authors & Date: Riccardo Rigon, February, 1999

References:  S. Wolfram, Mathematica 3.0, Cambridge University Press

Bugs & Limitations: Has been thought for functions not having vertical lines

*/
double interpolate(double x,DOUBLEMATRIX *cleandata,DOUBLEMATRIX* W);


/**
Name: interpolate_floatmatrix
 
Synopsis: interpolate_floatmatrix(FLOATMATRIX *matrice, float dt, float istante, FLOATVECTOR *vettore);

Description: linearly interpolate at the instant "istante" the values of a temporal series spaced by "dt".
		It is assumed that each row of a matrix of float contain a set of variables measured at the same time.

Inputs: 1) 	matrice: the data set that serve as basis for the interpolation; 
		2) 	dt: the temporal interval ;
		3) 	istante: the time when to interpolate the data;
		
Return:1) 	vettore: a vector with the interpolated values in the instant "istante".
			Its lenght equals the number of columns in the input matrix.

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Examples: APPLICATIONS/HYDROLOGY/geotop/geotop.c

See Also: interpolating_function

Authors & Date: Giacomo Bertoldi, September, 2000

References: 

Bugs & Limitations:
*/
void interpolate_floatmatrix(FLOATMATRIX *matrice, float dt, float istante, FLOATVECTOR *vettore);

/**
Name: interpolate_doublematrix
 
Synopsis: interpolate_doublematrix(DOUBLEMATRIX *matrice, float dt, float istante, DOUBLEVECTOR *vettore);

Description: linearly interpolate at the instant "istante" the values of a temporal series spaced by "dt".
		It is assumed that each row of a matrix of double contain a set of variables measured at the same time.

Inputs: 1) 	matrice: the data set that serve as basis for the interpolation; 
		2) 	dt: the temporal interval ;
		3) 	istante: the time when to interpolate the data;
		
Return:1) 	vettore: a vector with the interpolated values in the instant "istante".
			Its lenght equals the number of columns in the input matrix.

FILE: LIBRARIES/BASICMATHSTAT/t_datamanipulation.h, LIBRARIES/BASICMATHSTAT/datamanipulation.c

Examples: APPLICATIONS/HYDROLOGY/geotop/geotop.c

See Also: interpolating_function

Authors & Date: Giacomo Bertoldi, September, 2000

References: 

Bugs & Limitations:
*/
void interpolate_doublematrix(DOUBLEMATRIX *matrice, float dt, float istante, DOUBLEVECTOR *vettore);

/**


Name:   variance_doublematrix_column

Synopsis:  
double variance_doublematrix_column(DOUBLEMATRIX* net,long column,double mean);


Version: 0.96

Description: evaluate the variance of the set of data contained in the specified column 
of the specified doublematrix. If the field mean is set to 0, it returns the second moment


Authors &  Date: Riccardo Rigon, June 2000


Inputs: 1) the pointer to the matrix to be analyzed; 2) the number of the column 
to be analyzed; 3) the mean value used to evaluate the variance


Return: 

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

See Also: mean_doublematrix_column


*/

double variance_doublematrix_column(DOUBLEMATRIX* net,long column,double mean);


/**

Name:   mean_doublematrix_column

Synopsis:  
double mean_doublematrix_column(DOUBLEMATRIX* net,long column);


Version: 0.96

Description: evaluate the mean of the set of data contained in the specified column 
of the specified doublematrix. 


Authors &  Date: Riccardo Rigon, June 2000


Inputs: 1) the pointer to the matrix to be analyzed; 2) the number of the column 
to be analyzed; 


Return: 

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

See Also: variance_doublematrix_column

*/

double mean_doublematrix_column(DOUBLEMATRIX* net,long column);

/**

Name:   approximate_2_multiple

Synopsis:  
double approximate_2_multiple(double number,double div);

Version: 0.96

Description: approximate the 'number' to the closest smaller multiple of  'div'

Authors &  Date: Riccardo Rigon, January 2001


Inputs: 1) the number to be approximated; 2) the  smallest unit; 


Return: desired approximation

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

double approximate_2_multiple(double number,double div);

/**

Name:   ricampiona

Synopsis:  
DOUBLEMATRIX *ricampiona(DOUBLEMATRIX *cleandata, float ti2, float dt2, long n2);

Version: 0.96

Description: Ricampiona una serie di dati da dt inferiori a dt superiori

Authors &  Date: Giacomo Bertoldi, April 2001


Inputs: 1) the number to be approximated; 2) the  smallest unit; 


Return: desired approximation

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

DOUBLEMATRIX *ricampiona(DOUBLEMATRIX *cleandata, float ti2, float dt2, long n2, float dt1);




/*-------- quickinterpolate -----------------------------------------------

Synopsis:
  
DOUBLEMATRIX *quickinterpolate(DOUBLEMATRIX *cleandata, short intervals);
Version: 0.96

Description: Interpolate a matrix, giving a matrix with less rows

Authors &  Date: Giacomo Bertoldi, April 2001


Inputs: 1)	cleandata: matrix to interpolate 
		2)	intervals: number of intervals to aggregate


Return: interpolated matrix

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

DOUBLEMATRIX *quickinterpolate(DOUBLEMATRIX *cleandata, short intervals);


/**

Name:   mean_function

Synopsis:  
double mean_function(DOUBLEMATRIX *data, long n);


Version: 0.96

Description: evaluate the mean of the set of data contained in  the specified doublematrix:
the first column contain x - data,
the second column contain y - data.
The integral is evaluated with the trapezoidal rule.

Authors &  Date: Giacomo Bertoldi, April 2001


Inputs: 1) the pointer to the matrix to be analyzed;
		1) the number of rows of the matrix


Return: 

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

See Also: variance_doublematrix_column

*/

double mean_function(DOUBLEMATRIX *data, long n);

/**

Name:   cleandata_matrix

Synopsis:  
DOUBLEMATRIX *cleandata_matrix(DOUBLEMATRIX *cleandata, FLOATVECTOR *V, SHORTMATRIX *control);

Version: 0.96

Description: Elimina i novalues in una serie di dati interpolando linearmente tra i valori validi

Authors &  Date: Giacomo Bertoldi, April 2001


Inputs: 1) cleandata: the matrix to clean; 2) V: the vector of novalues; 
Outputs: 1) control: a matrix of the same size of cleandata with 1 for modified data, 0 for original data

Return: cleaned matrix

Bugs: not are allowed novalues in the first and the last row

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

DOUBLEMATRIX *cleandata_matrix(DOUBLEMATRIX *cleandata, FLOATVECTOR *V, SHORTMATRIX *control);

/*-------- fill_data ------------------------------------------------------------------*/

/**

Name:   fill_data

Synopsis:  
void fill_data(DOUBLEMATRIX *cleandata,long j,long n_i,long n_f,double x_i,double x_f);

Version: 0.96

Description: Elimina i novalues in una serie di dati interpolando linearmente tra i valori validi

Authors &  Date: Giacomo Bertoldi, April 2001


Inputs: 1) cleandata: the matrix to clean; 
	1) j: the colunm to clean
	2) n_i: the row before the first bad value
	3) n_l: the row past the last bad value
	4) x_i: the true value up
	5) x_i: the true value down
Outputs: 1) cleandata: the matrix to clean;
		2) control: a matrix of the same size of cleandata with 1 for modified data, 0 for original data


Bugs: not are allowed novalues in the first and the last row

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

void fill_data(DOUBLEMATRIX *cleandata,SHORTMATRIX *control,long j,long n_i,long n_f,double x_i,double x_f);

/*-------- aggregate ------------------------------------------------------------------*/

/**

Name:   aggregate

Synopsis:  
DOUBLEMATRIX *aggregate(DOUBLEMATRIX *data, long col)

Version: 0.96

Description: calculates the mean aggregating all data wich have the same value in the column col

Authors &  Date: Giacomo Bertoldi, October 2003


 input: data: double matrix with original data
		  col: the column with 
 return: double matrix aggregated data as mean values; 
   		   in the last column you have the number of aggregated element


Bugs: not are allowed novalues in the first and the last row

FIle: LIBRARIES/BASICSMATHSTAT/datamanipulation.c

*/

DOUBLEMATRIX *aggregate(DOUBLEMATRIX *data, long col, float nv);

