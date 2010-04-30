/**

Name: expdev

Synopsis: double expdev(long *idum);

Version: 1.0

Description: produces a random number with gaussian distribution


Authors & Date: Riccardo Rigon, 1998


Inputs: a pointer to a long  (if the long is negative the random number generator is
reinitialized)

Return: a random number with gaussian distribution (zero mean and standard deviation 1)

Notes: this code is copyrighted and you need NR to use it

Keywords: random numbers

References: Numerical Recipes, Second Edition, pag.287

*/
double expdev(long *idum);

/**

Name: poisdev

Synopsis:  double poisdev(float xm, long *idum);

Version: 1.0

Description: Generates a random number with poisson distribution


Inputs:  1- the distribution mean 2 - a pointer to a long  (if the long is negative the random number generator is
reinitialized)

Return: a randonm number with poisson distribution

See Also: expdev

Notes: this code is copyrighted and you need NR to use it

References: Numerical Recipes, Second Ediction, pag.294

*/

double poisdev(float xm, long *idum);

/**

Name: bnldev

Synopsis:  float bnldev(float pp, int n, long *idum);

Version: 1.0

Description: returns as a floating point number an integer value that is a random deviate drawn from a binomial
distribution of n trials each of probability pp

Inputs:  1- the probability of a trial 2 - the numbe of trials 3 - a pointer to a long  (if the long is negative the random number generator is
reinitialized)

Return: a randonm number with the above distribution

See Also: expdev

Notes: this code is copyrighted and you need NR to use it

References: Numerical Recipes, Second Ediction

*/
float bnldev(float pp, int n, long *idum);
/**

Name: gammln

Synopsis:  float gammln(float x);

Version: 1.0

Description: returns the value of  log(gamma(x))

Inputs:  1- the abscissa for which you need the value of the gamma

Return: a randonm number with the above distribution


Notes: this code is copyrighted and you need NR to use it

References: Numerical Recipes, Second Ediction, page 214

*/
float gammln(float xx);


double meandoublem(DOUBLEMATRIX *m,FLOATVECTOR *V);
double vardoublem(DOUBLEMATRIX *m,FLOATVECTOR *V);
