short transport_model(double Dt, ALLDATA *adt);
double Dt_transport_model(double Courant, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt);
void calculate_new_concentrations(double Dt, DOUBLEVECTOR *M, DOUBLEVECTOR *C, DOUBLEVECTOR *H, DOUBLEVECTOR *K, DOUBLEVECTOR *D, ALLDATA *adt);
void C_conversion_matrices_to_vector(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt);
void C_conversion_vector_to_matrices(DOUBLEVECTOR *V, DOUBLEMATRIX *Mland, DOUBLEMATRIX *Mchannel, ALLDATA *adt);

