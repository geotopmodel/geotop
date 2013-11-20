#include "constants.h"
#include "struct.geotop.h"
#include "pedo.funct.h"
#include "sparse_matrix.h"
#include "util_math.h"
#include "transport.model.h"
#include "meteodata.h"
#include <time.h>

extern long Nl, Nr, Nc; //max l, max r, max c
extern T_INIT *UV;	//dem information

//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************
//*****************************************************************************************************


short transport_model(double Dt, ALLDATA *adt){
/*
	if (Dt>Dt_transport_model(Courant, M, C, H, K, D, adt){
		return 1;
	}else{
		calculate_new_concentrations(Dt, M, C, H, K, D, adt);
		return 0;
	}
*/

//printf("\n\nHallo ich bin das Transport Modell\n\n");

return 0;

}

double Dt_transport_model(double Courant, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt){
	return 1.; //comment
}

void calculate_new_concentrations(double Dt, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt){}
void C_conversion_matrices_to_vector(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt){}
void C_conversion_vector_to_matrices(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt){}

/*


 ^-^
(0 0)
  v
(   )
 V V
 m m


*/

