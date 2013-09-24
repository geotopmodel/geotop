#include "rw_maps.h"

using namespace std;


//SHORTMATRIX *copyshort_doublematrix(DOUBLEMATRIX *M){
void copyshort_doublematrix(GeoMatrix<short>& S, GeoMatrix<double>& M){

//	SHORTMATRIX *S;
//	GeoMatrix<short> S;

	long r, c;

//	S=new_shortmatrix(M->nrh,M->nch);
	S.resize(M.getRows() ,M.getCols());
//	for(r=1;r<=M->nrh;r++){
	for(r=1;r<M.getRows();r++){
	//	for(c=1;c<=M->nch;c++){
		for(c=1;c<M.getCols();c++){
		//	S->co[r][c]=(short)M->co[r][c];
			S[r][c]=(short)M[r][c];
		}
	}

//	return(S);

}


//LONGMATRIX *copylong_doublematrix(DOUBLEMATRIX *M){
void copylong_doublematrix(GeoMatrix<long>& L, GeoMatrix<double>& M){

	long r, c;

//	L=new_longmatrix(M->nrh,M->nch);
	L.resize(M.getRows(),M.getCols());
//	for(r=1;r<=M->nrh;r++){
	for(r=1;r< M.getRows();r++){
	//	for(c=1;c<=M->nch;c++){
		for(c=1;c<M.getCols();c++){
		//	L->co[r][c]=(long)M->co[r][c];
			L[r][c]=(long)M[r][c];
		}
	}

}


//DOUBLEMATRIX *copydoublematrix_const(double c0, GeoMatrix<double>& Mref, double NOVALUE){
void copydoublematrix_const(double c0, GeoMatrix<double>& Mref,GeoMatrix<double>& M, double NOVALUE){

//	DOUBLEMATRIX *M;
	long r, c;

//	M=new_doublematrix(Mref->nrh,Mref->nch);
	M.resize(Mref.getRows(),Mref.getCols());
//	for(r=1;r<=M->nrh;r++){
	for(r=1;r<M.getRows();r++){
//		for(c=1;c<=M->nch;c++){
		for(c=1;c<M.getCols();c++){
		//	if(Mref->co[r][c]==NOVALUE){
			if(Mref[r][c]==NOVALUE){
			//	M->co[r][c]=NOVALUE;
				M[r][c]=NOVALUE;
			}else{
			//	M->co[r][c]=c0;
				M[r][c]=c0;
			}
		}
	}

//	return(M);

}
//----------------------

void write_suffix(char *suffix, long i, short start){

	short m=0, c=0, d=0, u=0;
	
	if(i>=0 && i<=9){
		m=0;
		c=0;
		d=0;
		u=i;
	}else if(i<=99){
		m=0;
		c=0;
		d=(short)(i/10.0);
		u=i-10.0*d;
	}else if(i<=999){
		m=0;
		c=(short)(i/100.0);
		d=(short)((i-100.0*c)/10.0);
		u=i-100.0*c-10*d;
	}else if(i<=9999){
		m=(short)(i/1000.0);
		c=(short)((i-1000.0*m)/100.0);
		d=(short)((i-1000.0*m-100.0*c)/10.0);
		u=i-1000*m-100.0*c-10*d;
	}else{
		t_error("Number too high");
	}
	
	m+=48;
	c+=48;
	d+=48;
	u+=48;

	suffix[start]=m;
	suffix[start+1]=c;
	suffix[start+2]=d;
	suffix[start+3]=u;
	
}
//overloaded function
void write_suffix(string suffix, long i, short start){

	short m=0, c=0, d=0, u=0;

	if(i>=0 && i<=9){
		m=0;
		c=0;
		d=0;
		u=i;
	}else if(i<=99){
		m=0;
		c=0;
		d=(short)(i/10.0);
		u=i-10.0*d;
	}else if(i<=999){
		m=0;
		c=(short)(i/100.0);
		d=(short)((i-100.0*c)/10.0);
		u=i-100.0*c-10*d;
	}else if(i<=9999){
		m=(short)(i/1000.0);
		c=(short)((i-1000.0*m)/100.0);
		d=(short)((i-1000.0*m-100.0*c)/10.0);
		u=i-1000*m-100.0*c-10*d;
	}else{
		t_error("Number too high");
	}

	m+=48;
	c+=48;
	d+=48;
	u+=48;

	suffix[start]=m;
	suffix[start+1]=c;
	suffix[start+2]=d;
	suffix[start+3]=u;

}

char *namefile_i(char *name, long i){

	char SSSS[ ]={"SSSS"};
	char *name_out;
	char *temp;
	
	write_suffix(SSSS, i, 0);
	
	temp=join_strings(name,SSSS);
	name_out=join_strings(temp,textfile);
	free(temp);
	
	return(name_out);
	
}


char *namefile_i_we(char *name, long i){

	char SSSS[ ]={"LSSSS"};
	char *name_out;
	
	write_suffix(SSSS, i, 1);
	
	name_out=join_strings(name,SSSS);
		
	return(name_out);
	
}	


char *namefile_i_we2(char *name, long i){

	char SSSS[ ]={"SSSS"};
	char *name_out;
	
	write_suffix(SSSS, i, 0);
	
	name_out=join_strings(name,SSSS);
		
	return(name_out);
	
}	

// TODO: Noori - supposed to return a pointer
//DOUBLEVECTOR *read_map_vector(short type, char *namefile, DOUBLEMATRIX *mask, T_INIT *grid, double no_value, LONGMATRIX *rc){
  GeoVector<double> read_map_vector(short type, char *namefile, GeoMatrix<double>& mask, TInit *grid, double no_value, GeoMatrix<long>& rc){
	
//	DOUBLEMATRIX *M;
	GeoMatrix<double> M;
//	DOUBLEVECTOR *V;
	GeoVector<double> V;
//	long i, n=rc->nrh;
	long i, n=rc.getRows()-1;
	
//	M = read_map(type, namefile, mask, grid, no_value);
	meteoio_readMap(string(namefile), M);

//	V = new_doublevector(n);
	V.resize(n);
	
	for (i=1; i<n; i++) {
	//	V->co[i] = M->co[rc->co[i][1]][rc->co[i][2]];
		V[i] = M[rc[i][1]][rc[i][2]];
	}
	
//	free_doublematrix(M);
	
	return V;
	
}



//void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV, long novalue){
void write_map(char *filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii(filename, type, M, UV, novalue);
		t_error("The Grass Ascii format not supported any more");
	}else if(format==3){
		write_esriascii(filename, type, M, UV, novalue);
	}

}
// overloaded function
void write_map(string filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii(filename, type, M, UV, novalue);
		t_error("The Grass Ascii format not supported any more");
	}else if(format==3){
		write_esriascii(filename, type, M, UV, novalue);
	}

}



//void write_map(char *filename, short type, short format, LONGMATRIX *M, T_INIT *UV, long novalue){
void write_map(char *filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii(filename, type, M, UV, novalue);
		t_error("The Grass Ascii format not supported any more");
	}else if(format==3){
		write_esriascii(filename, type, M, UV, novalue);
	}

}

//overloaded function
void write_map(string filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii(filename, type, M, UV, novalue);
		t_error("The Grass Ascii format not supported any more");
	}else if(format==3){
		write_esriascii(filename, type, M, UV, novalue);
	}

}

//void write_map_vector(char *filename, short type, short format, DOUBLEVECTOR *V, T_INIT *UV, long novalue, long **j, long nr, long nc){
/*
void write_map_vector(char *filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc){
	//	type=0  floating point
	//	type=1  integer

	//	format=1 fluidturtle
	//	format=2 grassascii
	//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii_vector(filename, type, V, j, nr, nc, UV, novalue);
		t_error("Grass ascii map format not supported any more");
	}else if(format==3){
		write_esriascii_vector(filename, type, V, j, nr, nc, UV, novalue);
	}

}
*/

//void write_map_vector(char *filename, short type, short format, DOUBLEVECTOR *V, T_INIT *UV, long novalue, long **j, long nr, long nc){
void write_map_vector(string filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc){
	//	type=0  floating point
	//	type=1  integer

	//	format=1 fluidturtle
	//	format=2 grassascii
	//	format=3 esriascii

	if(format==1){
		t_error("The fluidturtle format is not support any more");
	}else if(format==2){
		//write_grassascii_vector(filename, type, V, j, nr, nc, UV, novalue);
		t_error("Grass ascii map format not supported any more");
	}else if(format==3){
		write_esriascii_vector(filename, type, V, j, nr, nc, UV, novalue);
	}

}

//------------------------------------
//void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue){
void write_tensorseries(short a, long l, long i, char *filename, short type, short format, GeoTensor<double>& T, TInit *UV, long novalue){

//	a=0 non include "l" nel suffisso
//	a=1 include "l" nel suffisso
//	l:layer
//	i:temporal step
	
	char SSSSLLLLL[ ]={"SSSSLLLLL"};
	char SSSS[ ]={"SSSS"};		
	char *name;
	long r, c;
//	DOUBLEMATRIX *M;
	GeoMatrix<double> M;
		
	if(a==0){
		write_suffix(SSSS, i, 0);	
		name=join_strings(filename,SSSS);				
	}else if(a==1){
		write_suffix(SSSSLLLLL, i, 0);	
		write_suffix(SSSSLLLLL, l, 5);	
		name=join_strings(filename,SSSSLLLLL);		
	}else {
		t_error("Value not admitted");
	}

 
//	M=new_doublematrix(T->nrh,T->nch);
	M.resize(T.getRh(),T.getCh());
//	for(r=1;r<=T->nrh;r++){
	for(r=1;r<T.getRh();r++){
	//	for(c=1;c<=T->nch;c++){
		for(c=1;c<T.getCh();c++){
		//	M->co[r][c]=T->co[l][r][c];
			M[r][c]=T[l][r][c];
		}
	}

	write_map(name, type, format, M, UV, novalue);
	
//	free_doublematrix(M);
	free(name);

}


//void write_tensorseries_vector(short a, long l, long i, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc) {
void write_tensorseries_vector(short a, long l, long i, string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc){
	
	//	a=0 non include "l" nel suffisso
	//	a=1 include "l" nel suffisso
	//	l:layer
	//	i:temporal step
	
	char SSSSLLLLL[ ]={"SSSSLLLLL"};
	char SSSS[ ]={"SSSS"};		
//	char *name;
	string name;
//	long j, npoints=T->nch;
	long j, npoints=T.getCols();
//	DOUBLEVECTOR *V;
	GeoVector<double> V;
	
	if(a==0){
		write_suffix(SSSS, i, 0);	
	//	name=join_strings(filename,SSSS);
		name= filename+ SSSS;
	}else if(a==1){
		write_suffix(SSSSLLLLL, i, 0);	
		write_suffix(SSSSLLLLL, l, 5);	
	//	name=join_strings(filename,SSSSLLLLL);
		name= filename + SSSSLLLLL;
	}else {
		t_error("Value not admitted");
	}
	
	//V=new_doublevector(npoints);
	V.resize(npoints+1);

	for(j=1;j<=npoints-1;j++){
	//	V->co[j]=T->co[l][j];
		V[j]=T[l][j];
	}
	
	write_map_vector(name, type, format, V, UV, novalue, J, nr, nc);
	
	//free_doublevector(V);
//	free(name);
	
}


//---------------------------------------------

void write_tensorseries2(char *suf, long l, char *filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue){
	
	char LLLLL[ ]={"LLLLL"};
	char *temp1, *temp2;
	long r, c;
//	DOUBLEMATRIX *M;
	GeoMatrix<double> M;
		
	temp1 = join_strings(LLLLL, suf);
	write_suffix(temp1, l, 1);	
 
//	M = new_doublematrix(T->nrh,T->nch);
	M.resize(T->nrh+1,T->nch+1);

	for(r=1; r<=T->nrh; r++){
		for(c=1; c<=T->nch; c++){
		//	M->co[r][c] = T->co[l][r][c];
			M[r][c] = T->co[l][r][c];
		}
	}
	
	temp2 = join_strings(filename, temp1);
	write_map(temp2, type, format, M, UV, novalue);
	
//	free_doublematrix(M);
	free(temp1);
	free(temp2);

}

//void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc){
void write_tensorseries2_vector(char *suf, long l, char *filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc){
	
	char LLLLL[ ]={"LLLLL"};
//	char  *temp1,*temp2;
	string temp1, temp2;
//	long i, npoints=T->nch;
	long i, npoints=T.getCols();
//	DOUBLEVECTOR *V;
	GeoVector<double> V;
	
	temp1 = join_strings(LLLLL, suf);
	write_suffix(temp1, l, 1);	
	
	//V = new_doublevector(npoints);
	V.resize(npoints+1);
	
	for(i=1; i<=npoints; i++){
	//	V->co[i] = T->co[l][i];
		V[i] = T[l][i];
	}
	
//	temp2 = join_strings(filename, temp1);
	temp2 = filename + temp1;
	write_map_vector(temp2, type, format, V, UV, novalue, J, nr, nc);
	
	//free_doublevector(V);
//	free(temp1);
//	free(temp2);
	
}

//overloaded function
void write_tensorseries2_vector(string suf, long l, string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc){

	char LLLLL[ ]={"LLLLL"};
//	char *temp1, *temp2;
	string temp1, temp2;
//	long i, npoints=T->nch;
	long i, npoints=T.getCols();
//	DOUBLEVECTOR *V;
	GeoVector<double> V;

//	temp1 = join_strings(LLLLL, suf);
	temp1 = LLLLL + suf;
	write_suffix(temp1, l, 1);

	//V = new_doublevector(npoints);
	V.resize(npoints+1);

	for(i=1; i<npoints; i++){
	//	V->co[i] = T->co[l][i];
		V[i] = T[l][i];
	}

//	temp2 = join_strings(filename, temp1);
	temp2 = filename + temp1;
	write_map_vector(temp2, type, format, V, UV, novalue, J, nr, nc);

	//free_doublevector(V);
//	free(temp1);
//	free(temp2);

}


//---------------
//void write_tensorseries3(char *suffix, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV, long novalue){
void write_tensorseries3(char *suffix, char *filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue){

	long l;
	for(l=T->ndl;l<=T->ndh;l++){
		write_tensorseries2(suffix, l, filename, type, format, T, UV, novalue);
	}
}

//--------------------------
//void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, DOUBLEMATRIX *T, T_INIT *UV, long novalue, long **J, long nr, long nc){
void write_tensorseries3_vector(char *suffix, char *filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc){

	long l;
	// TODO: need check nrl
//	for(l=T->nrl;l<=T->nrh;l++){
	for(l=1;l<=T.getRows();l++){
		write_tensorseries2_vector(suffix, l, filename, type, format, T, UV, novalue, J, nr, nc);
	}
}
//overloaded function
void write_tensorseries3_vector(string suffix, string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc){

	long l;
	// TODO: need check nrl
//	for(l=T->nrl;l<=T->nrh;l++){
	for(l=1;l<T.getRows();l++){
		write_tensorseries2_vector(suffix, l, filename, type, format, T, UV, novalue, J, nr, nc);
	}
}



/*===================================================================================*/
/*===================functions copied from the file write_ascii.c====================*/
/*===================================================================================*/

//void write_esriascii(char *name, short type, LONGMATRIX *DTM, T_INIT *UV, long novalue){
void write_esriascii(char *name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;
	char *temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

	temp = join_strings(name,ascii_esri);
	f=fopen(temp,"w");

//	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"ncols         %ld\n",DTM.getCols());
//	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"nrows         %ld\n",DTM.getRows());
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

//	for(r=1;r<=DTM->nrh;r++){
	for(r=1;r<=DTM.getRows();r++){
	//	for(c=1;c<=DTM->nch;c++){
		for(c=1;c<=DTM.getCols();c++){
		//	if((long)DTM->co[r][c]==novalue){
			if((long)DTM[r][c]==novalue){
				fprintf(f,"%ld.0",novalue);
			}else{
				if(type==1){
				//	fprintf(f,"%ld",(long)(DTM->co[r][c]));
					fprintf(f,"%ld",(long)(DTM[r][c]));
				}else{
				//	fprintf(f,"%f",DTM->co[r][c]);
					fprintf(f,"%f",DTM[r][c]);
				}
			}
		//	if(c<DTM->nch) fprintf(f," ");
			if(c<DTM.getCols()) fprintf(f," ");
		}
	//	if(r<DTM->nrh) fprintf(f,"\n");
		if(r<DTM.getRows()) fprintf(f,"\n");
	}
	fprintf(f,"\n");// added by Matteo to avoid warnings when reading with R
	fclose(f);
	free(temp);
}
//overloaded function
void write_esriascii(string name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;
//	char *temp;
	string temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

//	temp = join_strings(name,ascii_esri);
	temp = name+ascii_esri;
	f=fopen(temp.c_str(),"w");

//	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"ncols         %ld\n",DTM.getCols());
//	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"nrows         %ld\n",DTM.getRows());
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

//	for(r=1;r<=DTM->nrh;r++){
	for(r=1;r<=DTM.getRows()-1;r++){
	//	for(c=1;c<=DTM->nch;c++){
		for(c=1;c<=DTM.getCols()-1;c++){
		//	if((long)DTM->co[r][c]==novalue){
			if((long)DTM[r][c]==novalue){
				fprintf(f,"%ld.0",novalue);
			}else{
				if(type==1){
				//	fprintf(f,"%ld",(long)(DTM->co[r][c]));
					fprintf(f,"%ld",(long)(DTM[r][c]));
				}else{
				//	fprintf(f,"%f",DTM->co[r][c]);
					fprintf(f,"%f",DTM[r][c]);
				}
			}
		//	if(c<DTM->nch) fprintf(f," ");
			if(c<DTM.getCols()) fprintf(f," ");
		}
	//	if(r<DTM->nrh) fprintf(f,"\n");
		if(r<DTM.getRows()) fprintf(f,"\n");
	}
	fprintf(f,"\n");// added by Matteo to avoid warnings when reading with R
	fclose(f);
//	free(temp);
}



//void write_esriascii(char *name, short type, DOUBLEMATRIX *DTM, T_INIT *UV, long novalue){
void write_esriascii(char *name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;
	char *temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

	temp = join_strings(name,ascii_esri);
	f=fopen(temp,"w");

//	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"ncols         %ld\n",DTM.getCols()-1);
//	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"nrows         %ld\n",DTM.getRows()-1);
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

//	for(r=1;r<=DTM->nrh;r++){
	for(r=1;r<DTM.getRows();r++){
	//	for(c=1;c<=DTM->nch;c++){
		for(c=1;c<DTM.getCols();c++){
		//	if((long)DTM->co[r][c]==novalue){
			if((long)DTM[r][c]==novalue){
				fprintf(f,"%ld.000",novalue);
			}else{
				if(type==1){
				//	fprintf(f,"%ld",(long)(DTM->co[r][c]));
					fprintf(f,"%ld",(long)(DTM[r][c]));
				}else{
				//	fprintf(f,"%f",DTM->co[r][c]);
					fprintf(f,"%.3f",DTM[r][c]);
				}
			}
		//	if(c<DTM->nch) fprintf(f," ");
			if(c<DTM.getCols()) fprintf(f," ");
		}
	//	if(r<DTM->nrh) fprintf(f,"\n");
		if(r<DTM.getRows()) fprintf(f,"\n");
	}
	fprintf(f,"\n");// added by Matteo to avoid warnings when reading with R
	fclose(f);
	free(temp);
}

//overloaded function
void write_esriascii(string name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue){

//	type=0  floating point
//	type=1  integer

	FILE *f;
	long r,c;
//	char *temp;
	string temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

//	temp = join_strings(name,ascii_esri);
	temp = name + ascii_esri;
	f=fopen(temp.c_str(),"w");

//	fprintf(f,"ncols         %ld\n",DTM->nch);
	fprintf(f,"ncols         %ld\n",DTM.getCols()-1);
//	fprintf(f,"nrows         %ld\n",DTM->nrh);
	fprintf(f,"nrows         %ld\n",DTM.getRows()-1);
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

//	for(r=1;r<=DTM->nrh;r++){
	for(r=1;r<DTM.getRows();r++){
	//	for(c=1;c<=DTM->nch;c++){
		for(c=1;c<DTM.getCols();c++){
		//	if((long)DTM->co[r][c]==novalue){
			if((long)DTM[r][c]==novalue){
				fprintf(f,"%ld.0",novalue);
			}else{
				if(type==1){
				//	fprintf(f,"%ld",(long)(DTM->co[r][c]));
					fprintf(f,"%ld",(long)(DTM[r][c]));
				}else{
				//	fprintf(f,"%f",DTM->co[r][c]);
					fprintf(f,"%.3f",DTM[r][c]);   // %.3f <|--- the 3 is how many decimal places to print.
				}
			}
		//	if(c<DTM->nch) fprintf(f," ");
			if(c<DTM.getCols()) fprintf(f," ");
		}
	//	if(r<DTM->nrh) fprintf(f,"\n");
		if(r<DTM.getRows()) fprintf(f,"\n");
	}
	fprintf(f,"\n");// added by Matteo to avoid warnings when reading with R
	fclose(f);
//	free(temp);
}


//===============

//void write_esriascii_vector(char *name, short type, DOUBLEVECTOR *DTM, long **j, long nr, long nc, T_INIT *UV, long novalue){
void write_esriascii_vector(char *name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue){
	//	type=0  floating point
	//	type=1  integer

	FILE *f;
	long r,c;
	char *temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

	temp = join_strings(name,ascii_esri);
	f=fopen(temp,"w");

	fprintf(f,"ncols         %ld\n",nc);
	fprintf(f,"nrows         %ld\n",nr);
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if (j[r][c] > 0) {
				if(type==1){
					//fprintf(f,"%ld",(long)(DTM->co[j[r][c]]));
					fprintf(f,"%ld",(long)(DTM[j[r][c]]));
				}else{
					//fprintf(f,"%f",DTM->co[j[r][c]]);
					fprintf(f,"%.3f",DTM[j[r][c]]);
				}
			}else {
				fprintf(f,"%ld.000",novalue);
			}
			if(c<nc) fprintf(f," ");
		}
		if(r<nr) fprintf(f,"\n");
	}
	fclose(f);
	free(temp);
}

void write_esriascii_vector(string name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue){
	//	type=0  floating point
	//	type=1  integer

    char *basedir;
    int ret = 0;

	FILE *f;
	long r,c;
//	char *temp;
    string temp;

//	if(UV->U->co[1]!=UV->U->co[2]){
	if(UV->U[1]!=UV->U[2]){
	//	printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U->co[2],UV->U->co[1]);
		printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n",UV->U[2],UV->U[1]);
		t_error("Fatal error");
	}

//	temp = join_strings(name ,ascii_esri);
	temp = name +ascii_esri;

    basedir = dirname(strdup(temp.c_str()));
    ret = mkdirp(basedir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(-1 == ret){
        t_error("write_esriascii_vector(): Unable to create parent directories of file" + temp);
    }
	f=fopen(temp.c_str(),"w");
    if(NULL == f){
        t_error("write_esriascii_vector(): Unable to open file `" + temp + "` in write mode.");
    }

	fprintf(f,"ncols         %ld\n",nc);
	fprintf(f,"nrows         %ld\n",nr);
	fprintf(f,"xllcorner     %f\n",UV->U[4]);
	fprintf(f,"yllcorner     %f\n",UV->U[3]);
	fprintf(f,"cellsize      %f\n",UV->U[1]);
	fprintf(f,"NODATA_value  %ld.0\n",novalue);

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if (j[r][c] > 0) {
				if(type==1){
				//	fprintf(f,"%ld",(long)(DTM->co[j[r][c]]));
					fprintf(f,"%ld",(long)(DTM[j[r][c]]));
				}else{
				//	fprintf(f,"%f",DTM->co[j[r][c]]);
					fprintf(f,"%.3f",DTM[j[r][c]]);
				}
			}else {

				fprintf(f,"%ld.000",novalue);
			}
			if(c<nc) fprintf(f," ");
		}
		if(r<nr) fprintf(f,"\n");
	}
	fclose(f);


//	free(temp);
}


//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

void error_message(short format, long n, long n1, long n2, long n3, char *name)
//format=1 grassascii
//format=2 esriascii

{
	if(n==n1 || n==n2 || n==n3){
		if(format==1) printf("File %s incompleted, end of file or end of line reached",join_strings(name,ascii_grass));
		if(format==2) printf("File %s incompleted, end of file or end of line reached",join_strings(name,ascii_esri));
		t_error("Fatal error");
	}
}


/*===================functions copied from geomorphology.0875.c============*/

//void curvature(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *c1, DOUBLEMATRIX *c2, DOUBLEMATRIX *c3, DOUBLEMATRIX *c4, long undef){
void curvature(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& c1, GeoMatrix<double>& c2, GeoMatrix<double>& c3, GeoMatrix<double>& c4, long undef){

	long r,c;
	long R1, R2, C1, C2;
//	long nc=topo->nch;
	long nc=topo.getCols()-1;
//	long nr=topo->nrh;
	long nr=topo.getRows()-1;
	double delta;

	// Compute the curvature.
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
		//	if((long)topo->co[r][c]!=undef){
			if((long)topo[r][c]!=undef){

			//	c1->co[r][c]=0.0;
				c1[r][c]=0.0;
				R1=r-1;
				R2=r+1;
				C1=c;
				C2=c;
				delta=deltay;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	c1->co[r][c] += (topo->co[R1][C1]+topo->co[R2][C2]-2.*topo->co[r][c])/pow(delta,2.);
						c1[r][c] += (topo[R1][C1]+topo[R2][C2]-2.*topo[r][c])/pow(delta,2.);
					}
				}

			//	c2->co[r][c]=0.0;
				c2[r][c]=0.0;
				R1=r;
				R2=r;
				C1=c+1;
				C2=c-1;
				delta=deltax;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	c2->co[r][c] += (topo->co[R1][C1]+topo->co[R2][C2]-2.*topo->co[r][c])/pow(delta,2.);
						c2[r][c] += (topo[R1][C1]+topo[R2][C2]-2.*topo[r][c])/pow(delta,2.);
					}
				}

			//	c3->co[r][c]=0.0;
				c3[r][c]=0.0;
				R1=r-1;
				R2=r+1;
				C1=c-1;
				C2=c+1;
				delta=sqrt(deltax*deltay);
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	c3->co[r][c] += (topo->co[R1][C1]+topo->co[R2][C2]-2.*topo->co[r][c])/pow(delta,2.);
						c3[r][c] += (topo[R1][C1]+topo[R2][C2]-2.*topo[r][c])/pow(delta,2.);
					}
				}

			//	c4->co[r][c]=0.0;
				c4[r][c]=0.0;
				R1=r-1;
				R2=r+1;
				C1=c+1;
				C2=c-1;
				delta=sqrt(deltax*deltay);
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	c4->co[r][c] += (topo->co[R1][C1]+topo->co[R2][C2]-2.*topo->co[r][c])/pow(delta,2.);
						c4[r][c] += (topo[R1][C1]+topo[R2][C2]-2.*topo[r][c])/pow(delta,2.);
					}
				}


			}else {

			//	c1->co[r][c] = (double)undef;
				c1[r][c] = (double)undef;
			//	c2->co[r][c] = (double)undef;
				c2[r][c] = (double)undef;
			//	c3->co[r][c] = (double)undef;
				c3[r][c] = (double)undef;
			//	c4->co[r][c] = (double)undef;
				c4[r][c] = (double)undef;

			}

		}
	}
}

//------------------------
//short is_boundary(long r, long c, DOUBLEMATRIX *dem, long novalue){
short is_boundary(long r, long c, GeoMatrix<double>& dem, long novalue){

	long ir, ic;
	short yes = 0;

	ir=-1;
	ic=0;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=-1;
	ic=1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=0;
	ic=1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=1;
	ic=1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=1;
	ic=0;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=1;
	ic=-1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=0;
	ic=-1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	ir=-1;
	ic=-1;
//	if( (long)dem->co[r+ir][c+ic]==novalue ) yes = 1;
	if( (long)dem[r+ir][c+ic]==novalue ) yes = 1;

	return yes;

}


//long row(double N, long nrows, T_INIT *UV, long novalue){
long row(double N, long nrows, TInit *UV, long novalue){

	long cnt;

//	if(N<UV->U->co[3] || N>UV->U->co[3]+nrows*UV->U->co[1]){
	if(N<UV->U[3] || N>UV->U[3]+nrows*UV->U[1]){
		return novalue;
	}else {
		cnt=0;
		do{
			cnt++;
	//	}while(UV->U->co[3]+(nrows-cnt)*UV->U->co[1]>N);
		}while(UV->U[3]+(nrows-cnt)*UV->U[1]>N);
		return cnt;
	}
}


//long col(double E, long ncols, T_INIT *UV, long novalue){
long col(double E, long ncols, TInit *UV, long novalue){

	long cnt;

//	if(E<UV->U->co[4] || E>UV->U->co[4]+ncols*UV->U->co[2]){
	if(E<UV->U[4] || E>UV->U[4]+ncols*UV->U[2]){
		return novalue;
	}else{
		cnt=0;
		do{
			cnt++;
	//	}while(UV->U->co[4]+cnt*UV->U->co[2]<E);
		}while(UV->U[4]+cnt*UV->U[2]<E);
		return cnt;
	}
}



//Presa da geomorphology099 e modificato der_min
//void nablaquadro_mask(DOUBLEMATRIX *Z0,SHORTMATRIX *curv,DOUBLEVECTOR *U,DOUBLEVECTOR *V)
void nablaquadro_mask(GeoMatrix<double>& Z0, GeoMatrix<short>& curv, GeoVector<double>& U, GeoVector<double>& V)

{

short y;
long i,j,h,rows,cols;
double grid[9],z[9],derivate2;
double der_min=0.00001; /*limite per la limite per la planarita'*/

short v[13][2] = {           { 0, 0},
                             { 0, 1},
                             {-1, 1},
                             {-1, 0},
                             {-1,-1},
                             { 0,-1},
                             { 1,-1},
                             { 1, 0},
                             { 1, 1},
                             { 0, 0},
                             { 0, 0},
                             { 0, 0},
                             { 0, 0}             };

grid[0]=0;
//grid[1]=grid[5]=U->co[1];
grid[1]=grid[5]=U[1];
//grid[3]=grid[7]=U->co[2];
grid[3]=grid[7]=U[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

//rows=Z0->nrh;
rows=Z0.getRows()-1;
//cols=Z0->nch;
cols=Z0.getCols()-1;

for(i=2;i<=rows-1;i++){
   	for(j=2;j<=cols-1;j++){
   	//	z[0]=Z0->co[i][j];
   		z[0]=Z0[i][j];
   	//	if(z[0]!=V->co[2]){
   		if(z[0]!=V[2]){
       		y=1;
        	for(h=1;h<=8;h++){
        	//	z[h]=Z0->co[i+v[h][0]][j+v[h][1]];
        		z[h]=Z0[i+v[h][0]][j+v[h][1]];
			//	if(z[h]==V->co[2]){
        		if(z[h]==V[2]){
					y=0;
		        	break;
        		}
       		}
        	if(y==0){
        	//	curv->co[i][j]=1;
        		curv[i][j]=1;
        	}else{
		    	derivate2=0.5*((z[1]+z[5]-2*z[0])/(grid[1]*grid[1])+ (z[3]+z[7]-2*z[0])/(grid[3]*grid[3]));
		    	derivate2=derivate2+0.5*((z[2]+z[4]+z[6]+z[8]-4*z[0])/(grid[6]*grid[6]));

		    	if(fabs(derivate2)<=der_min || derivate2>der_min){  //plane or concave
		    	//	curv->co[i][j]=0;
		    		curv[i][j]=0;
		    	}else{
		    	//	curv->co[i][j]=1;	//convex
		    		curv[i][j]=1;		//convex
		    	}
		 	}
    	}
	}
}
}

/*====================copied function from init.c ==================*/

//void initmatrix(double val, DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue){
void initmatrix(double val, GeoMatrix<double>& destination, GeoMatrix<double>& origin, double novalue){

	long r,c;
//	for(r=1;r<=destination->nrh;r++){
	for(r=1;r<destination.getRows();r++){
	//	for(c=1;c<=destination->nch;c++){
		for(c=1;c<destination.getCols();c++){
		//	if(origin->co[r][c]!=novalue) destination->co[r][c]=val;
			if(origin[r][c]!=novalue) destination[r][c]=val;
		}
	}
}



//----------------------------------
