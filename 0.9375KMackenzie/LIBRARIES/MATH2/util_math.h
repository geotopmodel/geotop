/*----------------------------------------------------------------------------------------------------------*/
/*Used to solve tridiagonal linear system*/
/*----------------------------------------------------------------------------------------------------------*/


				
void tridiag(short a, long r, long c, long nx, DOUBLEVECTOR *diag_inf, DOUBLEVECTOR *diag, DOUBLEVECTOR *diag_sup, DOUBLEVECTOR *b, DOUBLEVECTOR *e);

double norm(DOUBLEVECTOR *V);
