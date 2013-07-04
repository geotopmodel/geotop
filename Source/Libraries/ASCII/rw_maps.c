
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
    
#include "turtle.h"
#include "rw_maps.h"
#include "import_ascii.h"
#include "write_ascii.h"
#include "extensions.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "tensor3D.h"


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

SHORTMATRIX *copyshort_doublematrix(DOUBLEMATRIX *M){

	SHORTMATRIX *S;
	long r, c;
	
	S=new_shortmatrix(M->nrh,M->nch);
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			S->co[r][c]=(short)M->co[r][c];
		}
	}
	
	return(S);

}
	
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

LONGMATRIX *copylong_doublematrix(DOUBLEMATRIX *M){

	LONGMATRIX *L;
	long r, c;
	
	L=new_longmatrix(M->nrh,M->nch);
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			L->co[r][c]=(long)M->co[r][c];
		}
	}
	
	return(L);

}
	
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *copydouble_shortmatrix(SHORTMATRIX *S){

	DOUBLEMATRIX *M;
	long r, c;
	
	M=new_doublematrix(S->nrh,S->nch);
	for(r=1;r<=S->nrh;r++){
		for(c=1;c<=S->nch;c++){
			M->co[r][c]=(double)S->co[r][c];
		}
	}
	
	return(M);

}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *copydouble_longmatrix(LONGMATRIX *L){

	DOUBLEMATRIX *M;
	long r, c;
	
	M=new_doublematrix(L->nrh,L->nch);
	for(r=1;r<=L->nrh;r++){
		for(c=1;c<=L->nch;c++){
			M->co[r][c]=(double)L->co[r][c];
		}
	}
	
	return(M);

}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *copydoublematrix_const(double c0, DOUBLEMATRIX *Mref, double NOVALUE){

	DOUBLEMATRIX *M;
	long r, c;
	
	M=new_doublematrix(Mref->nrh,Mref->nch);
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(Mref->co[r][c]==NOVALUE){
				M->co[r][c]=NOVALUE;
			}else{
				M->co[r][c]=c0;
			}
		}
	}
	
	return(M);
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *multiplydoublematrix(double f, DOUBLEMATRIX *Mref, double NOVALUE){

	DOUBLEMATRIX *M;
	long r, c;
	
	M=new_doublematrix(Mref->nrh,Mref->nch);
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			if(Mref->co[r][c]==NOVALUE){
				M->co[r][c]=NOVALUE;
			}else{
				M->co[r][c]=f*Mref->co[r][c];
			}
		}
	}
	
	return(M);
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void build_doubletensor(DOUBLETENSOR *T, DOUBLEMATRIX *M, long l){

	long r,c;
	
	if(l<=0 || l>T->ndh) t_error("Invalid doubletensor construction");
	if(T->nrh!=M->nrh) t_error("Invalid doubletensor construction");
	if(T->nch!=M->nch) t_error("Invalid doubletensor construction");	
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			T->co[l][r][c]=M->co[r][c];
		}
	}
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *extract_doublematrix(DOUBLETENSOR *T, long l){

	long r,c;
	DOUBLEMATRIX *M;
	
	if(l<=0 || l>T->ndh) t_error("Invalid doubletensor extraction");
	
	M=new_doublematrix(T->nrh,T->nch);
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			M->co[r][c]=T->co[l][r][c];
		}
	}
	
	return(M);
	
}



//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *extract_fromtensor(DOUBLETENSOR *T, long l){

	long r, c;
	DOUBLEMATRIX *M;
	
	if(l<1 || l>T->ndh) t_error("cannot extract a matrix from a tensor");
	
	M=new_doublematrix(T->nrh,T->nch);
	for(r=1;r<=T->nrh;r++){
		for(c=1;c<=T->nch;c++){
			M->co[r][c]=T->co[l][r][c];
		}
	}
	
	return(M);
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLETENSOR *build_frommatrix(DOUBLEMATRIX *M, long l, long lmax){

	long r,c;
	DOUBLETENSOR *T;
	
	if(lmax<1) t_error("cannot build a tensor from a matrix");
	if(l>lmax) t_error("cannot build a tensor from a matrix");
	
	T=new_doubletensor(lmax, M->nrh, M->nch);
	initialize_doubletensor(T, 0.0);
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			T->co[l][r][c]=M->co[r][c];
		}
	}	
	
	return(T);
	
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_frommatrix(long l, DOUBLEMATRIX *M, DOUBLETENSOR *T){

	long r,c;
	
	if(l<0 || l>T->ndh) t_error("cannot write a tensor from a matrix");
	if(M->nrh!=T->nrh) t_error("cannot write a tensor from a matrix");
	if(M->nch!=T->nch) t_error("cannot write a tensor from a matrix");
	
	for(r=1;r<=M->nrh;r++){
		for(c=1;c<=M->nch;c++){
			T->co[l][r][c]=M->co[r][c];
		}
	}		
	
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void fmultiplydoublematrix(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double f, double novalue){
	
	long r,c;
	
	if(origin->nrh!=destination->nrh) t_error("cannot copy matrix");
	if(origin->nch!=destination->nch) t_error("cannot copy matrix");
	for(r=1;r<=origin->nrh;r++){
		for(c=1;c<=origin->nch;c++){
			if(origin->co[r][c]!=novalue){
				destination->co[r][c]=f*origin->co[r][c];
			}else{
				destination->co[r][c]=novalue;
			}
		}
	}
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void assignnovalue(DOUBLEMATRIX *destination, DOUBLEMATRIX *origin, double novalue){
	
	long r,c;
	
	if(origin->nrh!=destination->nrh) t_error("cannot assign novalue");
	if(origin->nch!=destination->nch) t_error("cannot assign novalue");
	for(r=1;r<=origin->nrh;r++){
		for(c=1;c<=origin->nch;c++){
			if(origin->co[r][c]==novalue) destination->co[r][c]=novalue;
		}
	}
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_suffix(char *suffix, long i, short start){

	short m, c, d, u;
	
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

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

short existing_file(char *name){

	//if the file exists gives 1 (fluidturtle), 2(grassascii), 3(esriascii), 0 if the file doesn't exist

	short a=0;
	FILE *f;
	char *ft,*esri,*grass;
	
	ft=join_strings(name,ascii_ft);
	esri=join_strings(name,ascii_esri);
	grass=join_strings(name,ascii_grass);
	
	if( (f=fopen(ft,"r"))!=NULL ){
		a=1;
		fclose(f);
	}else if( (f=fopen(grass,"r"))!=NULL ){
		a=2;
		fclose(f);	
	}else if( (f=fopen(esri,"r"))!=NULL ){
		a=3;
		fclose(f);
	}
	
	free(ft);
	free(esri);
	free(grass);
	
	return(a);

}

short existing_file_text(char *name){

	//if the file exists gives 1, 0 if the file doesn't exist

	short a=0;
	FILE *f;
	char *temp;
	
	temp=join_strings(name,textfile);

	if( (f=fopen(temp,"r"))!=NULL ){
		a=1;
		fclose(f);	
	}
	
	free(temp);
	
	return(a);

}
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *read_map(short a, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref){

//	a=0 non usa Mref, UVref output
//	a=1 non esegue controllo non values, Mref e UVref input
//	a=2 esegue controllo novalues, Mref e UVref input

	DOUBLEMATRIX *M;
	FILE *f;
	long index, r, c, nr, nc;
	double *header, *m;
	double Dxmap, Dymap, X0map, Y0map, novalue;
	char *temp;
	T_INIT *UV;
		
	if (existing_file(filename)==1){
		
		if(a>0){
			UV=(T_INIT *)malloc(sizeof(T_INIT));
			if(!UV) t_error("UV was not allocated in read_map");
		}

		temp=join_strings(filename,ascii_ft);		
		f=t_fopen(temp,"r");
		index=read_index(f,PRINT);
		
		if(a==0){
			UVref->U=read_doublearray(f,PRINT);
			UVref->V=read_doublearray(f,PRINT);	
		}else{
			UV->U=read_doublearray(f,PRINT);
			UV->V=read_doublearray(f,PRINT);
		}
		
		M=read_doublematrix(f,"a",PRINT);
		t_fclose(f);
				
		if(a>0){
			//Check header
			if(UV->U->co[2]!=UVref->U->co[2]){
				printf("Dx in %s file is different from Dx in DTM file! \n",temp);
				t_error("Inconsistent map");
			}
			if(UV->U->co[1]!=UVref->U->co[1]){
				printf("Dy in %s file is different from Dy in DTM file! \n",temp);
				t_error("Inconsistent map");
			}		
			if(UV->U->co[4]!=UVref->U->co[4]){
				printf("X0 in %s file is different from X0 in DTM file! \n",temp);
				t_error("Inconsistent map");
			}		
			if(UV->U->co[3]!=UVref->U->co[3]){
				printf("Y0 in %s file is different from Y0 in DTM file! \n",temp);
				t_error("Inconsistent map");
			}		
			if(M->nrh!=Mref->nrh){
				printf("Number of rows in %s file is not consistent with DTM file! \n",temp);
				t_error("Inconsistent map");
			}
			if(M->nch!=Mref->nch){
				printf("Number of columns in %s file is not consistent with DTM file! \n",temp);
				t_error("Inconsistent map");
			}
		
			//The novalue is imposed egual to that of DTM file:*/
			for(r=1;r<=M->nrh;r++){
				for(c=1;c<=M->nch;c++){
					if(M->co[r][c]==UV->V->co[2]) M->co[r][c]=UVref->V->co[2];
					if(a>1){
						if(M->co[r][c]==UVref->V->co[2] && Mref->co[r][c]!=UVref->V->co[2]){
							printf("Novalues not consistent in %s file (1)",temp);
							printf("\nr:%ld c:%ld Mref:%f M:%f -> set at 0 \n",r,c,Mref->co[r][c], M->co[r][c]);	
							M->co[r][c] = 0.0;						
							//t_error("Inconsistent map");
						}
						if(M->co[r][c]!=UVref->V->co[2] && Mref->co[r][c]==UVref->V->co[2]) M->co[r][c]=UVref->V->co[2];
					}else{
						if(M->co[r][c]==UVref->V->co[2] && Mref->co[r][c]!=UVref->V->co[2]){
							printf("Novalues not consistent in %s file",temp);
							printf("\nr:%ld c:%ld Mref:%f M:%f\n",r,c,Mref->co[r][c], M->co[r][c]);
							stop_execution();
						}
					}
				}
			}
		
			free_doublevector(UV->U);
			free_doublevector(UV->V);
			free(UV);
			
		}

	}else{	
		if(a>0){
			novalue=UVref->V->co[2];
		}else{
			novalue=-9999.0;
		}
		header = (double *) malloc(6*sizeof(double));
						
		if(existing_file(filename)==2){	
			m=read_grassascii(header, novalue, filename);
			nr=(long)header[4];
			nc=(long)header[5];		
			Dxmap=(header[2]-header[3])/((long)header[5]);
			Dymap=(header[0]-header[1])/((long)header[4]);
			X0map=header[3];
			Y0map=header[1];
			temp=join_strings(filename,ascii_grass);
		}else if(existing_file(filename)==3){
			m=read_esriascii(header, novalue, filename);
			nr=(long)header[1];
			nc=(long)header[0];				
			Dxmap=header[4];
			Dymap=header[4];
			X0map=header[2];
			Y0map=header[3];
			temp=join_strings(filename,ascii_esri);
		}else{
			printf("The file %s doesn't exist\n",filename);
			t_error("Fatal error");
		}
		
		free(header);
		if(a>0){
			//Check header	
			if(Dxmap!=UVref->U->co[2]){
				printf("Dx in %s file is different from Dx in DTM file! \n",temp);
				t_error("Inconsistent map");
			}
			if(Dymap!=UVref->U->co[1]){
				printf("Dy:%f in %s file is different from Dy:%f in DTM file! \n",Dymap,temp,UVref->U->co[1]);
				t_error("Inconsistent map");
			}
			if(X0map!=UVref->U->co[4]){
				printf("X0 in %s file is different from X0 in DTM file! \n",temp);
				t_error("Inconsistent map");
			}
			if(Y0map!=UVref->U->co[3]){	
				printf("Y0 in %s file is different from Y0 in DTM file! \n",temp);
				t_error("Inconsistent map");
			}
			if(nr!=Mref->nrh){
				printf("Number of rows in %s file (%ld) is not consistent with DTM file (%ld)! \n",temp,nr,Mref->nrh);
				t_error("Inconsistent map");
			}
			if(nc!=Mref->nch){
				printf("Number of columns in %s file is not consistent with DTM file! \n",temp);
				t_error("Inconsistent map");
			}

		}else{
			UVref->U=new_doublevector(4);
			UVref->V=new_doublevector(2);
			UVref->U->co[2]=Dxmap;
			UVref->U->co[1]=Dymap;
			UVref->U->co[4]=X0map;
			UVref->U->co[3]=Y0map;
			UVref->V->co[1]=-1;
			UVref->V->co[2]=novalue;
		}		
				
		//assign values and check novalues
		M=new_doublematrix(nr,nc);
		for(r=1;r<=nr;r++){
			for(c=1;c<=nc;c++){
				M->co[r][c]=m[(r-1)*nc+c-1];
				if(a>1){
					if (M->co[r][c]==UVref->V->co[2] && Mref->co[r][c]!=UVref->V->co[2]){
						printf("Novalues not consistent in %s file (2)",temp);
						printf("\nr:%ld c:%ld Mref:%f M:%f -> set at 0 \n",r,c,Mref->co[r][c], M->co[r][c]);	
						M->co[r][c] = 0.0;						
						//t_error("Inconsistent map");	
					}
					if(M->co[r][c]!=UVref->V->co[2] && Mref->co[r][c]==UVref->V->co[2]) M->co[r][c]=UVref->V->co[2];
				}
			}
		}
		
		free(m);
	}
	
	free(temp);
	
	return(M);
	
}
	
	
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLEMATRIX *read_mapseries(long i, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref){

	char SSSS[ ]={"LSSSS"};
	char *name;
	DOUBLEMATRIX *M;

	write_suffix(SSSS, i, 1);	
	name=join_strings(filename,SSSS);
	M=read_map(2, name, Mref, UVref);	
	free(name);

	return(M);
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

DOUBLETENSOR *read_tensor(long nl, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref){
	long l;	
	DOUBLEMATRIX *M;
	DOUBLETENSOR *T;
			
	if(nl<1) t_error("The dimension of the tensor must be greater than or equal to 1");
	
	for(l=1;l<=nl;l++){
		M=read_mapseries(l, filename, Mref, UVref);
		if(l==1){
			T=build_frommatrix(M, l, nl);
		}else{
			write_frommatrix(l, M, T);
		}
		free_doublematrix(M);
	}
	
	return(T);
}



//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------		

DOUBLETENSOR *read_maptensor(long i, long lmax, char *filename, DOUBLEMATRIX *Mref, T_INIT *UVref){

	char SSSSLLLLL[ ]={"SSSSLLLLL"};
	char *name;
	DOUBLETENSOR *T;
	DOUBLEMATRIX *M;
	long l;
	
	if(i<0 || lmax<1) t_error("cannot read a tensor with null or negative columns");

	write_suffix(SSSSLLLLL, i, 0);	
	
	for(l=1;l<=lmax;l++){
		write_suffix(SSSSLLLLL, l, 5);
		name=join_strings(filename,SSSSLLLLL);
		M=read_map(2, name, Mref, UVref);
		free(name);
		if(l==1){
			T=build_frommatrix(M, l, lmax);
		}else{
			write_frommatrix(l, M, T);
		}
	}

	return(T);
	
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_map(char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV){

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii
	
	if(format==1){
		write_fluidturtle(filename, type, M, UV);
	}else if(format==2){
		write_grassascii(filename, type, M, UV);
	}else if(format==3){
		write_esriascii(filename, type, M, UV);
	}
	
}	

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_mapseries(long i, char *filename, short type, short format, DOUBLEMATRIX *M, T_INIT *UV){
	
	char SSSS[ ]={"SSSS"};	
	char *name;
	
	write_suffix(SSSS, i, 0);	
	name=join_strings(filename, SSSS);
	write_map(name, type, format, M, UV);
	free(name);

}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_tensorseries(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV){

//	a=0 non include "l" nel suffisso
//	a=1 include "l" nel suffisso
//	l:layer
//	i:temporal step
	
	char SSSSLLLLL[ ]={"SSSSLLLLL"};
	char SSSS[ ]={"SSSS"};		
	char *name;
	long r, c;
	DOUBLEMATRIX *M;
		
	if(a==0){
		write_suffix(SSSS, i, 0);	
		name=join_strings(filename,SSSS);				
	}else if(a==1){
		write_suffix(SSSSLLLLL, i, 0);	
		write_suffix(SSSSLLLLL, l, 5);	
		name=join_strings(filename,SSSSLLLLL);		
	}
 
	M=new_doublematrix(T->nrh,T->nch);
	for(r=1;r<=T->nrh;r++){
		for(c=1;c<=T->nch;c++){
			M->co[r][c]=T->co[l][r][c];
		}
	}

	write_map(name, type, format, M, UV);
	
	free_doublematrix(M);
	free(name);

}

void write_tensorseries_bis(short a, long l, long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV){

//	a=0 non include "l" nel suffisso
//	a=1 include "l" nel suffisso
//	l:layer
//	i:temporal step
	
	char SSSSLLLLL[ ]={"LLLLLNNNNN"};
	char SSSS[ ]={"NNNN"};		
	char *name;
	long r, c;
	DOUBLEMATRIX *M;
		
	if(a==0){
		write_suffix(SSSS, i, 0);	
		name=join_strings(filename,SSSS);				
	}else if(a==1){
		write_suffix(SSSSLLLLL, l, 1);	
		write_suffix(SSSSLLLLL, i, 6);	
		name=join_strings(filename,SSSSLLLLL);		
	}
 
	M=new_doublematrix(T->nrh,T->nch);
	

	for(r=1;r<=T->nrh;r++){
		for(c=1;c<=T->nch;c++){
			M->co[r][c]=T->co[l][r][c];
		}
	}
	
	write_map(name, type, format, M, UV);
	
	free_doublematrix(M);
	free(name);

}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void write_tensorseries2(long i, char *filename, short type, short format, DOUBLETENSOR *T, T_INIT *UV){

	long l;
	for(l=T->ndl;l<=T->ndh;l++){
		write_tensorseries_bis(1, l, i, filename, type, format, T, UV);
	}
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

