#include <string>
#include "t_io.h"

#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <iostream>

/*WORKING_POSITION=SEEK_SET;*/

FILE * EXTERNAL_FILE;

LONGVECTOR *EXTERNAL_P;

char  EXTERNAL_FILE_NAME[256],OLD_NAME[256];

long int EXTERNAL_FILE_POSITION=SEEK_SET;

HEADER EXTERNAL_HEADER;



short OPENYES=0;

long IO_FILES_COUNTER=0;
long IO_STRINGS_COUNTER=0;
long  WORKING_POSITION=0;
long IO_PARMS_COUNTER=0;


t_keywords T_KEYWORDS={{"2","ascii","binary"},

							  {"7","char","short","int","long","float","double","string"},

							  {"5","array","vector","matrix","list","tensor"},

							  {"2","->","<-"},

							  {"2"," ",","},

							  {"2","{","}"}};


/*-----------------------------------------------------------------------*/

int mkdirp(const char *pathname, mode_t mode)
{
    /* From http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/ */
    char parent[PATH_MAX], *p;

    /* make a parent directory path */
    strncpy(parent, pathname, sizeof(parent));
    parent[sizeof(parent) - 1] = '\0';
    for(p = parent + strlen(parent); *p != '/' && p != parent; p--);
    *p = '\0';

    /* try make parent directory */
    if(p != parent && mkdirp(parent, mode) != 0)
        return -1;

    /* make this one if parent has been made */
    if(mkdir(pathname, mode) == 0)
        return 0;

    /* if it already exists that is fine */
    if(errno == EEXIST){
        struct stat sb;
        stat(pathname, &sb);
        if(!S_ISDIR(sb.st_mode)){
            /* pathname is NOT a directory! */
            return -1;
        }
        return 0;
    }
    return -1;
}

//void t_error(char *error_text)
void t_error(std::string error_text)
/* Error handling */
{

	/* void exit(); */
	fprintf(stderr,"\nError::Run time error\n");
	fprintf(stderr,"Error::%s\n",error_text.c_str());
	fprintf(stderr,"........exiting...\n");
	exit(1);
}


/**-----------------------------------------------------------------------*/
FILE *t_fopen(const char *name,const char *mode)
/* Safe file open */
{

//char newname[256];
FILE *fp=NULL;

char *lBasePtr = strdup(name) ;
char *basedir = dirname(lBasePtr);
free(lBasePtr) ;
int ret = 0;
    
ret = mkdirp(basedir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
if(-1 == ret){
    fprintf(stderr, "ERROR: Unable to create parent directory `%s`. Exiting.\n", basedir);
    exit(1);
}

if((fp=fopen(name,mode))==NULL){
	printf("%s",name);    
	t_error(" The specified file could not be opened ");

	return NULL;
	
}else{

	return fp;

}

}

/**-----------------------------------------------------------------------*/
FILE * t_fclose(FILE * stream)
/* Safe file close */
{


if(stream==NULL){
	printf(" An attemp was made to close an already closed file ");
}else{		
	 fclose(stream);
}

	 return NULL;

}

/*
TODO: removeme
char *get_workingdirectory(void )

{

    //char buffer[64*FILENAME_MAX];
    char *bf=NULL,*pathfile="$WorkingPath";
    //long len;
    long i;
    short a;
    FILE *istream;

    istream=fopen(pathfile,"r");
    if(istream){

        i=0;
        a=0;
        do{
            if(i==0){
                bf = (char *) malloc(sizeof(char));
            }else{
                bf = (char *)realloc(bf,(i+1)*sizeof(char));
            }
            bf[i]=fgetc(istream);
            if(bf[i]==10 || bf[i]==-1){
                a=1;
                bf[i]=0;
            }
            i+=1;
        }while(a==0);

        fclose(istream);

    }else{

//		printf("ENTER THE WORKING DIRECTORY PATH:\n");
//		scanf("%s",&buffer);
//		len=64*FILENAME_MAX;
//
//    	if(len > (64*FILENAME_MAX)){
//			t_error("Maximum path length exceeded");
//		} else {
//			bf=(char *)malloc(len*sizeof(char));
//		}
//
//		strcpy(bf,buffer);

        t_error("You have to specify aworking directory when you run the executable");
    }
    return bf;
}
*/

std::string get_workingdirectory()
{
    std::string retString ;
    char const * const pathfile = getenv("WorkingPath");

    if(pathfile == NULL)
    {
        t_error("You need to export the WorkingPath environment variable before to run the program");
        return retString ;
    }

    std::ifstream inputFile ;
    inputFile.open (pathfile, std::ifstream::in);
    std::getline (inputFile, retString);

    if( retString == "" )
    {
        t_error(std::string("no valid path found on file: ") + pathfile);
    }

    return retString;
}

/**-------------------------------------------------------------------------*/

/*
char *join_strings(char const * const first, char const * const second)

{

	char *string;
	int len=strlen(first)+strlen(second)+2;

	string=(char*)malloc(len*sizeof(char));
	string=strcpy(string,first);
	string=strcat(string,second);

	return string;
}
*/

/*================functions copied from datamanipulation.c===========*/

/*---------------------------------------------------------------------------*/
void initialize_longvector(LONGVECTOR *L, long sign)

{

	long i;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nh; i++) {
				L->co[i] = sign;
			}
		} else {
			t_error("This longvector was no properly allocated");
		}
	} else {
		t_error("A null vector was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_shortvector(SHORTVECTOR *L, short sign)

{

	long i;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nh; i++) {
				L->co[i] = sign;
			}
		} else {
			t_error("This shortvector was no properly allocated");
		}
	} else {
		t_error("A null vector was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_doublevector(DOUBLEVECTOR *L, double sign)

{

	long i;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = L->nl; i <= L->nh; i++) {
				L->co[i] = sign;
			}
		} else {
			t_error("This doublevector was no properly allocated");
		}
	} else {
		t_error("A null vector was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_floatvector(FLOATVECTOR *L, float sign)

{

	long i;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nh; i++) {
				L->co[i] = sign;
			}
		} else {
			t_error("This floatvector was no properly allocated");
		}
	} else {
		t_error("A null vector was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_shortmatrix(SHORTMATRIX *L, short sign)

{

	long i, j;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nrh; i++) {
				for (j = 1; j <= L->nch; j++) {
					L->co[i][j] = sign;
				}
			}
		} else {
			t_error("This matrix was no properly allocated");
		}
	} else {
		t_error("A null matrix was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_longmatrix(LONGMATRIX *L, long sign)

{

	long i, j;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nrh; i++) {
				for (j = 1; j <= L->nch; j++) {
					L->co[i][j] = sign;
				}
			}
		} else {
			t_error("This matrix was no properly allocated");
		}
	} else {
		t_error("A null matrix was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_floatmatrix(FLOATMATRIX *L, float sign)

{

	long i, j;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nrh; i++) {
				for (j = 1; j <= L->nch; j++) {
					L->co[i][j] = sign;
				}
			}
		} else {
			t_error("This matrix was no properly allocated");
		}
	} else {
		t_error("A null matrix was addressed");
	}
}

/*---------------------------------------------------------------------------*/
void initialize_doublematrix(DOUBLEMATRIX *L, double sign)

{

	long i, j;

	if (L != NULL) {
		if (L->isdynamic == 1) {
			for (i = 1; i <= L->nrh; i++) {
				for (j = 1; j <= L->nch; j++) {
					L->co[i][j] = sign;
				}
			}
		} else {
			t_error("This matrix was no properly allocated");
		}
	} else {
		t_error("A null matrix was addressed");
	}
}

/*--------------------------------------------------------------------------*/

void copy_shortmatrix(SHORTMATRIX *origin, SHORTMATRIX *destination)

{

	long i, j;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A matrix was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nrh < 1 || destination->nrh < 1 || origin->nch < 1
			|| destination->nch < 1) {

		t_error("A matrix was not allocated properly");

	} else if (origin->nrh != destination->nrh || origin->nch
			!= destination->nch) {

		t_error("The matrixes do not have the same dimensions");

	}

	for (i = 1; i <= origin->nrh; i++) {
		for (j = 1; j <= origin->nch; j++) {

			destination->co[i][j] = origin->co[i][j];

		}
	}

}

/*--------------------------------------------------------------------------*/

void copy_intmatrix(INTMATRIX *origin, INTMATRIX *destination)

{

	long i, j;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A matrix was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nrh < 1 || destination->nrh < 1 || origin->nch < 1
			|| destination->nch < 1) {

		t_error("A matrix was not allocated properly");

	} else if (origin->nrh != destination->nrh || origin->nch
			!= destination->nch) {

		t_error("The matrixes do not have the same dimensions");

	}

	for (i = 1; i <= origin->nrh; i++) {
		for (j = 1; j <= origin->nch; j++) {

			destination->co[i][j] = origin->co[i][j];

		}
	}

}

/*--------------------------------------------------------------------------*/

void copy_longmatrix(LONGMATRIX *origin, LONGMATRIX *destination)

{

	long i, j;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A matrix was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nrh < 1 || destination->nrh < 1 || origin->nch < 1
			|| destination->nch < 1) {

		t_error("A matrix was not allocated properly");

	} else if (origin->nrh != destination->nrh || origin->nch
			!= destination->nch) {

		t_error("The matrixes do not have the same dimensions");

	}

	for (i = 1; i <= origin->nrh; i++) {
		for (j = 1; j <= origin->nch; j++) {

			destination->co[i][j] = origin->co[i][j];

		}
	}

}

/*--------------------------------------------------------------------------*/

void copy_floatmatrix(FLOATMATRIX *origin, FLOATMATRIX *destination)

{

	long i, j;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A matrix was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nrh < 1 || destination->nrh < 1 || origin->nch < 1
			|| destination->nch < 1) {

		t_error("A matrix was not allocated properly");

	} else if (origin->nrh != destination->nrh || origin->nch
			!= destination->nch) {

		t_error("The matrixes do not have the same dimensions");

	}

	for (i = 1; i <= origin->nrh; i++) {
		for (j = 1; j <= origin->nch; j++) {

			destination->co[i][j] = origin->co[i][j];

		}
	}

}

/*--------------------------------------------------------------------------*/

void copy_doublematrix(DOUBLEMATRIX *origin, DOUBLEMATRIX *destination)

{

	long i, j;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A matrix was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nrh < 1 || destination->nrh < 1 || origin->nch < 1
			|| destination->nch < 1) {

		t_error("A matrix was not allocated properly");

	} else if (origin->nrh != destination->nrh || origin->nch
			!= destination->nch) {

		t_error("The matrixes do not have the same dimensions");

	}

	for (i = 1; i <= origin->nrh; i++) {
		for (j = 1; j <= origin->nch; j++) {

			destination->co[i][j] = origin->co[i][j];

		}
	}

}

/*--------------------------------------------------------------------------*/

void copy_shortvector(SHORTVECTOR *origin, SHORTVECTOR *destination)

{

	long i;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A vector was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nh < 1 || destination->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (origin->nh != destination->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= origin->nh; i++) {

		destination->co[i] = origin->co[i];

	}

}

/*--------------------------------------------------------------------------*/

void copy_intvector(INTVECTOR *origin, INTVECTOR *destination)

{

	long i;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A vector was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nh < 1 || destination->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (origin->nh != destination->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= origin->nh; i++) {

		destination->co[i] = origin->co[i];

	}

}

/*--------------------------------------------------------------------------*/

void copy_longvector(LONGVECTOR *origin, LONGVECTOR *destination)

{

	long i;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A vector was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nh < 1 || destination->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (origin->nh != destination->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= origin->nh; i++) {

		destination->co[i] = origin->co[i];

	}

}

/*--------------------------------------------------------------------------*/

void copy_floatvector(FLOATVECTOR *origin, FLOATVECTOR *destination)

{

	long i;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A vector was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nh < 1 || destination->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (origin->nh != destination->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= origin->nh; i++) {

		destination->co[i] = origin->co[i];

	}

}

/*--------------------------------------------------------------------------*/

//void copy_doublevector(DOUBLEVECTOR *origin, DOUBLEVECTOR *destination)
void copy_doublevector(DOUBLEVECTOR *origin, DOUBLEVECTOR *destination)

{

	long i;

	if (origin == NULL || destination == NULL || origin->co == NULL
			|| destination->co == NULL) {

		t_error("A vector was not allocated");

	} else if (origin->isdynamic != 1 || destination->isdynamic != 1
			|| origin->nh < 1 || destination->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (origin->nh != destination->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= origin->nh; i++) {

		destination->co[i] = origin->co[i];

	}

}

/*--------------------------------------------------------------------------*/

void add_doublevector(DOUBLEVECTOR *small1, DOUBLEVECTOR *big1)

{

	long i;

	if (small1 == NULL || big1 == NULL || small1->co == NULL || big1->co == NULL) {

		t_error("A vector was not allocated");

	} else if (small1->isdynamic != 1 || big1->isdynamic != 1 || small1->nh < 1
			|| big1->nh < 1) {

		t_error("A vector was not allocated properly");

	} else if (small1->nh != big1->nh) {

		t_error("The vector do not have the same dimensions");

	}

	for (i = 1; i <= small1->nh; i++) {

		big1->co[i] += small1->co[i];

	}

}

void print_doublematrix_elements(DOUBLEMATRIX *m,long maxcols)
/* Write a matrix of  double to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%f ",m->co[i][j]);
		    if(j%maxcols==0 && j!= m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');
}
}
//----------------------------
//void multipass_topofilter(long ntimes, DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout, long novalue, long n){
void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n){

    long i ;
//	DOUBLEMATRIX *M;
	GeoMatrix<double> M;

//	M=new_doublematrix(Zin->nrh, Zin->nch);

////for (r=1; r<=Zin->nrh; r++) {
//	for (r=1; r<Zin.getRows(); r++) {
//	//	for (c=1; c<=Zin->nch; c++) {
//		for (c=1; c<Zin.getCols(); c++) {
//		//	Zout->co[r][c] = Zin->co[r][c];
//			Zout[r][c] = Zin[r][c];
//		}
//	}

	Zout=Zin;

	for (i=1; i<=ntimes; i++) {

//		for (r=1; r<=Zout->nrh; r++) {
//			for (c=1; c<=Zout->nch; c++) {
//				M->co[r][c] = Zout->co[r][c];
//			}
//		}
    M=Zout;

		topofilter(M, Zout, novalue, n);
	}
//	free_doublematrix(M);
}

//----------------------------------------
//void find_min_max(DOUBLEMATRIX *M, long novalue, double *max, double *min){
void find_min_max(GeoMatrix<double>& M, long novalue, double *max, double *min){

//	long r, c, nr=M->nrh, nc=M->nch;
	long r, c, nr=M.getRows(), nc=M.getCols();

	*max=-1.E99;
	*min=1.E99;

//	for (r=1; r<=nr; r++) {
	for (r=1; r<nr; r++) {
	//	for (c=1; c<=nc; c++) {
		for (c=1; c<nc; c++) {
		//	if ((long)M->co[r][c] != novalue) {
			if ((long)M[r][c] != novalue) {
			//	if(*max < M->co[r][c]) *max = M->co[r][c];
				if(*max < M[r][c]) *max = M[r][c];
			//	if(*min > M->co[r][c]) *min = M->co[r][c];
				if(*min > M[r][c]) *min = M[r][c];
			}
		}
	}
}

//-----------------------------------------
//double norm_1(DOUBLEVECTOR *V, long nbeg, long nend){
double norm_1(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;
	for(l=nbeg;l<=nend;l++){
	//	N += fabs(V->co[l]);
		N += fabs(V[l]);
	}

	return(N);

}

//***************************************************************************


/* Calculation of the sky view factor for each pixels:
   Input:  - N         number of parts in which one wants to divide the horizon
           - UV        format file with the dimension of pixel and the novalue
           - input     matrix with elevation (DTM)
           - convess   matrix with concave zones 0 and covex zones 1
   Output: - (sky)     matrix with sky view factor
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program sky of Pegoretti
                                                 */

//void sky_view_factor(DOUBLEMATRIX *sky, long N, T_INIT *UV, DOUBLEMATRIX *input, SHORTMATRIX *convess, long novalue)
void sky_view_factor(GeoMatrix<double>& sky, long N, TInit *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue)
{
 long i,j,t,m,n,p,q,h,k,r,s; //counter
 double deltateta; //amplitude of the angles in which the horizon is divided
 DOUBLEMATRIX *alfa; //matrices with the angles of the direction
 DOUBLEVECTOR *vv; //vector with the view factor of the current pixel for one of the N parts
 DOUBLEVECTOR *v; //vector with the minimum view factor of the current pixel for one of the N parts
 double vvv; //mean of the sky view for a pixel of the N parts
 long nr=input.getRows()-1;
 long nc=input.getCols()-1;

//	 if(sky->nrh!=input->nrh) t_error("Sky view factor fatal error, number of rows not consistent");
 	 if(sky.getRows()-1!=nr) t_error("Sky view factor fatal error, number of rows not consistent");
// 	 if(sky->nch!=input->nch) t_error("Sky view factor fatal error, number of cols not consistent");
 	 if(sky.getCols()-1!=nc) t_error("Sky view factor fatal error, number of cols not consistent");

// Computation of the matrix with the angles of the direction
// alfa=new_doublematrix(2*input->nrh-1,2*input->nch-1);
   alfa=new_doublematrix(2*nr-1,2*nc-1);
   initialize_doublematrix(alfa,(double)novalue); //initialisation with novalue

// for(i=1;i<=2*input->nrh-1;i++){
   for(i=1;i<=2*nr-1;i++){
   //  for(j=1;j<=2*input->nch-1;j++){
	   for(j=1;j<=2*nc-1;j++){
       //  if(i<=input->nrh && j<input->nch){
		   if(i<=nr && j<nc){
           //  alfa->co[i][j]=3.0/2.0*GTConst::Pi+atan(((input->nrh-i)*UV->U->co[1])/((input->nch-j)*UV->U->co[1]));
			   alfa->co[i][j]=3.0/2.0*GTConst::Pi+atan(((nr-i)*UV->U[1])/((nc-j)*UV->U[1]));
            }
         // if(i>input->nrh && j<=input->nch){
		    if(i>nr && j<=nc){
            //  alfa->co[i][j]=GTConst::Pi+atan(((input->nch-j)*UV->U->co[1])/((i-input->nrh)*UV->U->co[1]));
		    	 alfa->co[i][j]=GTConst::Pi+atan(((nc-j)*UV->U[1])/((i-nr)*UV->U[1]));
            }
         // if(i>=input->nrh && j>input->nch){
		    if(i>=nr && j>nc){
            //  alfa->co[i][j]=GTConst::Pi/2.0+atan(((i-input->nrh)*UV->U->co[1])/
            //                        ((j-input->nch)*UV->U->co[1]));
		        alfa->co[i][j]=GTConst::Pi/2.0+atan(((i-nr)*UV->U[1])/
		    	                                ((j-nc)*UV->U[1]));

            }
        //  if(i<input->nrh && j>=input->nch){
		    if(i<nr && j>=nc){
           //   alfa->co[i][j]=atan(((j-input->nch)*UV->U->co[1])/((input->nrh-i)*UV->U->co[1]));
		    	alfa->co[i][j]=atan(((j-nc)*UV->U[1])/((nr-i)*UV->U[1]));
            }
      }
 }

 // Computation of matrix with sky view factor:
// for(i=1;i<=sky->nrh;i++){
   for(i=1;i< sky.getRows();i++){
	  //for(j=1;j<=sky->nch;j++){
		for(j=1;j<sky.getCols();j++){
	//	sky->co[i][j]=(double)novalue;
		sky[i][j]=(double)novalue;
	}
 }

 v=new_doublevector(N);
 vv=new_doublevector(N);
 deltateta=2.0*GTConst::Pi/N;

// for(i=1;i<=input->nrh;i++){
   for(i=1;i<=nr;i++){
//	for(j=1;j<=input->nch;j++){
	for(j=1;j<=nc;j++){
	//	if ((long)input->co[i][j]!=novalue){ //computation only of novalue pixels
		if ((long)input[i][j]!=novalue){ //computation only of novalue pixels
		   	for(t=1;t<=N;t++){
		       	v->co[t]=1.0;
		   	}
		//  m=input->nrh-i+1;
		   	m=nr-i+1;
		// 	n=input->nch-j+1;
		 	n=nc-j+1;
		// 	p=m+input->nrh-1;
			p=m+nr-1;
		//  q=n+input->nch-1;
			q=n+nc-1;
           	for(h=m;h<=p;h++){
		      	for(k=n;k<=q;k++){
		         	for(t=1;t<=N;t++){
		               	if (alfa->co[h][k]>=(t-1)*deltateta && alfa->co[h][k]<t*deltateta){
		                  	r=h-m+1;
		                  	s=k-n+1;
		                 // if (convess->co[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
		                    if (convess[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
//		                        vv->co[t]=1-sin(atan((input->co[r][s]-input->co[i][j])
//		                         	          /(sqrt(pow((r-i),2)+pow((s-j),2))*UV->U->co[1])));
		                        vv->co[t]=1-sin(atan((input[r][s]-input[i][j])
		         		                         	          /(sqrt(pow((r-i),2)+pow((s-j),2))*UV->U[1])));
                             	if (vv->co[t]<v->co[t]){
                                   v->co[t]=vv->co[t];
                             	}
                          	}
                          	break;
                       	}
                 	}
              	}
           	}
           	vvv=0.0;
           	for(t=1;t<=N;t++){
               	vvv=vvv+v->co[t];
           	}
        //  sky->co[i][j]=(1.0/N*vvv);
           	sky[i][j]=(1.0/N*vvv);
         }
   	}
  // 	printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",100.0*(double)i/(double)sky->nrh);
 		printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",100.0*(double)i/(double)sky.getRows()-1);
}

free_doublematrix(alfa);
free_doublevector(v);
free_doublevector(vv);

}

//***************************************************************************
//DOUBLEMATRIX *find_aspect(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef){
void find_aspect(DOUBLEMATRIX *topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M){

	long r, c;
	long nc=topo->nch;
	long nr=topo->nrh;
//	DOUBLEMATRIX *M;
//	M=new_doublematrix(nr, nc);

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if((long)topo->co[r][c]!=undef){
			//	M->co[r][c] = (180./GTConst::Pi) * (3.0 / 2.0 * GTConst::Pi - atan2(dzdy->co[r][c], dzdx->co[r][c]));
				M[r][c] = (180./GTConst::Pi) * (3.0 / 2.0 * GTConst::Pi - atan2(dzdy[r][c], dzdx[r][c]));
				if (M[r][c]>=360.0) M[r][c] -= 360.0;
			}else {
				M[r][c] = (double)undef;
			}
		}
	}
}

//overloaded function noori
void find_aspect(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M){

	long r, c;
//	long nc=topo->nch;
	long nc=topo.getCols();
//	long nr=topo->nrh;
	long nr=topo.getRows();
//	DOUBLEMATRIX *M;
//	M=new_doublematrix(nr, nc);
	M.resize(nr,nc,-9999.0);
//	for(r=1;r<=nr;r++){
	for(r=1;r<nr;r++){
	//	for(c=1;c<=nc;c++){
		for(c=1;c<nc;c++){
		//	if((long)topo->co[r][c]!=undef){
			if((long)topo[r][c]!=undef){
			//	M->co[r][c] = (180./GTConst::Pi) * (3.0 / 2.0 * GTConst::Pi - atan2(dzdy->co[r][c], dzdx->co[r][c]));
				M[r][c] = (180./GTConst::Pi) * (3.0 / 2.0 * GTConst::Pi - atan2(dzdy[r][c], dzdx[r][c]));
				if (M[r][c]>=360.0) M[r][c] -= 360.0;
			}else {
				M[r][c] = (double)undef;
			}
		}
	}
}

//---------------
void topofilter(DOUBLEMATRIX *Zin, DOUBLEMATRIX *Zout, long novalue, long n){

	long r, c, nr, nc, ir, ic, i;
	DOUBLEVECTOR *values;
	long cnt;

	values=new_doublevector((2*n+1)*(2*n+1));

	nr=Zin->nrh;
	nc=Zin->nch;

	for (r=1; r<=nr; r++) {
		for (c=1; c<=nc; c++) {
			if ((long)Zin->co[r][c] != novalue) {
				cnt=0;
				for (ir=-n; ir<=n; ir++) {
					for (ic=-n; ic<=n; ic++) {
						if (r+ir>=1 && r+ir<=nr && c+ic>=1 && c+ic<=nc) {
							if((long)Zin->co[r+ir][c+ic]!=novalue){
								cnt++;
								values->co[cnt]=Zin->co[r+ir][c+ic];
							}
						}
					}
				}

				/*order_values(values, cnt);
				 if (fmod(cnt, 2)>1.E-5) {
				 Zout->co[r][c] = values->co[(long)(0.5*(cnt+1))];
				 }else{
				 Zout->co[r][c] = 0.5*(values->co[(long)(0.5*(cnt))]+values->co[(long)(0.5*(cnt)+1)]);
				 }*/

				Zout->co[r][c] = 0.;
				for (i=1; i<=cnt; i++) {
					Zout->co[r][c] += values->co[i]/(double)cnt;
				}

			}else {

				Zout->co[r][c] = (double)novalue;

			}

		}
	}

	free_doublevector(values);
}
// overloade function noori
void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n){

	long r, c, nr, nc, ir, ic, i;
	DOUBLEVECTOR *values;
	long cnt;

	values=new_doublevector((2*n+1)*(2*n+1));

//	nr=Zin->nrh;
	nr=Zin.getRows();
//	nc=Zin->nch;
	nc=Zin.getCols();

//	for (r=1; r<=nr; r++) {
	for (r=1; r<nr; r++) {
	//	for (c=1; c<=nc; c++) {
		for (c=1; c<nc; c++) {
		//	if ((long)Zin->co[r][c] != novalue) {
			if ((long)Zin[r][c] != novalue) {
				cnt=0;
				for (ir=-n; ir<=n; ir++) {
					for (ic=-n; ic<=n; ic++) {
						if (r+ir>=1 && r+ir<=nr && c+ic>=1 && c+ic<=nc) {
						//	if((long)Zin->co[r+ir][c+ic]!=novalue){
							if((long)Zin[r+ir][c+ic]!=novalue){
								cnt++;
							//	values->co[cnt]=Zin->co[r+ir][c+ic];
								values->co[cnt]=Zin[r+ir][c+ic];
							}
						}
					}
				}

				/*order_values(values, cnt);
				 if (fmod(cnt, 2)>1.E-5) {
				 Zout->co[r][c] = values->co[(long)(0.5*(cnt+1))];
				 }else{
				 Zout->co[r][c] = 0.5*(values->co[(long)(0.5*(cnt))]+values->co[(long)(0.5*(cnt)+1)]);
				 }*/

			//	Zout->co[r][c] = 0.;
				Zout[r][c] = 0.;
				for (i=1; i<=cnt; i++) {
				//	Zout->co[r][c] += values->co[i]/(double)cnt;
					Zout[r][c] += values->co[i]/(double)cnt;
				}

			}else {

			//	Zout->co[r][c] = (double)novalue;
				Zout[r][c] = (double)novalue;

			}

		}
	}

	free_doublevector(values);
}

//------------------------
void find_slope(double deltax, double deltay, DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef){

	long r,c,R1,C1,R2,C2;
	long nc=topo->nch;
	long nr=topo->nrh;
	double a, delta;

	// Find dzdx.
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if((long)topo->co[r][c]!=undef){

				a=0.0;
				R1=r;
				R2=r;
				C1=c-1;
				C2=c+1;
				delta=deltax;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
					if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
						if( (topo->co[R2][C2] - topo->co[r][c]) * (topo->co[r][c] - topo->co[R1][C1]) < 0){
							if( fabs(topo->co[r][c] - topo->co[R1][C1]) > fabs(topo->co[R2][C2] - topo->co[r][c]) ){
								a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
							}else {
								a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
							}
						}else {
							a += (topo->co[R2][C2] - topo->co[R1][C1]) / (2.*delta);
						}
					}else if((long)topo->co[R1][C1]==undef && (long)topo->co[R2][C2]!=undef){
						a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
					}else if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]==undef){
						a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
					}
				}
				dzdx->co[r][c] = a;

				a=0.0;
				R1=r+1;
				R2=r-1;
				C1=c;
				C2=c;
				delta=deltay;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
					if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
						if( (topo->co[R2][C2] - topo->co[r][c]) * (topo->co[r][c] - topo->co[R1][C1]) < 0){
							if( fabs(topo->co[r][c] - topo->co[R1][C1]) > fabs(topo->co[R2][C2] - topo->co[r][c]) ){
								a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
							}else {
								a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
							}
						}else {
							a += (topo->co[R2][C2] - topo->co[R1][C1]) / (2.*delta);
						}
					}else if((long)topo->co[R1][C1]==undef && (long)topo->co[R2][C2]!=undef){
						a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
					}else if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]==undef){
						a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
					}
				}
				dzdy->co[r][c] = a;

				//Some compilers will not allow dzdx and dzdy to both be 0.0 in the atan2 computation.
				//if (fabs(dzdx->co[r][c])<1.E-10) dzdx->co[r][c] = 1.E-10;
				if (fabs(dzdy->co[r][c])<1.E-10) dzdy->co[r][c] = 1.E-10;

			}else {

				dzdx->co[r][c] = (double)undef;
				dzdy->co[r][c] = (double)undef;

			}

		}
	}
}


//--------
// overloaded function noori
void find_slope(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef){

	long r,c,R1,C1,R2,C2;
//	long nc=topo->nch;
	long nc=topo.getCols() - 1;
//	long nr=topo->nrh;
	long nr=topo.getRows() - 1;
	double a, delta;

	// Find dzdx.
	for(r=1;r<nr;r++){
		for(c=1;c<nc;c++){
		//	if((long)topo->co[r][c]!=undef){
			if((long)topo[r][c]!=undef){

				a=0.0;
				R1=r;
				R2=r;
				C1=c-1;
				C2=c+1;
				delta=deltax;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	if( (topo->co[R2][C2] - topo->co[r][c]) * (topo->co[r][c] - topo->co[R1][C1]) < 0){
						if( (topo[R2][C2] - topo[r][c]) * (topo[r][c] - topo[R1][C1]) < 0){
						//	if( fabs(topo->co[r][c] - topo->co[R1][C1]) > fabs(topo->co[R2][C2] - topo->co[r][c]) ){
							if( fabs(topo[r][c] - topo[R1][C1]) > fabs(topo[R2][C2] - topo[r][c]) ){
							//	a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
								a += (topo[r][c] - topo[R1][C1]) / delta;
							}else {
							//	a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
								a += (topo[R2][C2] - topo[r][c]) / delta;
							}
						}else {
						//	a += (topo->co[R2][C2] - topo->co[R1][C1]) / (2.*delta);
							a += (topo[R2][C2] - topo[R1][C1]) / (2.*delta);
						}
				//	}else if((long)topo->co[R1][C1]==undef && (long)topo->co[R2][C2]!=undef){
					}else if((long)topo[R1][C1]==undef && (long)topo[R2][C2]!=undef){
					//	a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
						a += (topo[R2][C2] - topo[r][c]) / delta;
				//	}else if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]==undef){
					}else if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]==undef){
					//	a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
						a += (topo[r][c] - topo[R1][C1]) / delta;
					}
				}
			//	dzdx->co[r][c] = a;
				dzdx[r][c] = a;

				a=0.0;
				R1=r+1;
				R2=r-1;
				C1=c;
				C2=c;
				delta=deltay;
				if(R1>=1 && R1<=nr && R2>=1 && R2<=nr && C1>=1 && C1<=nc && C2>=1 && C2<=nc){
				//	if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]!=undef){
					if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]!=undef){
					//	if( (topo->co[R2][C2] - topo->co[r][c]) * (topo->co[r][c] - topo->co[R1][C1]) < 0){
						if( (topo[R2][C2] - topo[r][c]) * (topo[r][c] - topo[R1][C1]) < 0){
						//	if( fabs(topo->co[r][c] - topo->co[R1][C1]) > fabs(topo->co[R2][C2] - topo->co[r][c]) ){
							if( fabs(topo[r][c] - topo[R1][C1]) > fabs(topo[R2][C2] - topo[r][c]) ){
							//	a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
								a += (topo[r][c] - topo[R1][C1]) / delta;
							}else {
							//	a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
								a += (topo[R2][C2] - topo[r][c]) / delta;
							}
						}else {
						//	a += (topo->co[R2][C2] - topo->co[R1][C1]) / (2.*delta);
							a += (topo[R2][C2] - topo[R1][C1]) / (2.*delta);
						}
				//	}else if((long)topo->co[R1][C1]==undef && (long)topo->co[R2][C2]!=undef){
					}else if((long)topo[R1][C1]==undef && (long)topo[R2][C2]!=undef){
					//	a += (topo->co[R2][C2] - topo->co[r][c]) / delta;
						a += (topo[R2][C2] - topo[r][c]) / delta;
				//	}else if((long)topo->co[R1][C1]!=undef && (long)topo->co[R2][C2]==undef){
					}else if((long)topo[R1][C1]!=undef && (long)topo[R2][C2]==undef){
					//	a += (topo->co[r][c] - topo->co[R1][C1]) / delta;
						a += (topo[r][c] - topo[R1][C1]) / delta;
					}
				}
			//	dzdy->co[r][c] = a;
				dzdy[r][c] = a;

				//Some compilers will not allow dzdx and dzdy to both be 0.0 in the atan2 computation.
				//if (fabs(dzdx->co[r][c])<1.E-10) dzdx->co[r][c] = 1.E-10;
				if (fabs(dzdy[r][c])<1.E-10) dzdy[r][c] = 1.E-10;

			}else {

			//	dzdx->co[r][c] = (double)undef;
				dzdx[r][c] = (double)undef;
			//	dzdy->co[r][c] = (double)undef;
				dzdy[r][c] = (double)undef;

			}

		}
	}
}


//---------------------------------

DOUBLEMATRIX *find_max_slope(DOUBLEMATRIX *topo, DOUBLEMATRIX *dzdx, DOUBLEMATRIX *dzdy, long undef){

	long r, c;
	long nc=topo->nch;
	long nr=topo->nrh;
	DOUBLEMATRIX *M;

	M=new_doublematrix(nr, nc);

	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if((long)topo->co[r][c]!=undef){
				M->co[r][c] = (180./GTConst::Pi) * atan(pow(pow(dzdx->co[r][c], 2.0) + pow(dzdy->co[r][c], 2.0), 0.5));
			}else {
				M->co[r][c] = (double)undef;
			}

		}
	}

	return M;
}


//overloaded function noori
void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M){

	long r, c;
//	long nc=topo->nch;
	long nc=topo.getCols();
//	long nr=topo->nrh;
	long nr=topo.getRows();
//	DOUBLEMATRIX *M;

//	M=new_doublematrix(nr, nc);
	M.resize(nr, nc);

//	for(r=1;r<=nr;r++){
	for(r=1;r<nr;r++){
	//	for(c=1;c<=nc;c++){
		for(c=1;c<nc;c++){
		//	if((long)topo->co[r][c]!=undef){
			if((long)topo[r][c]!=undef){
			//	M->co[r][c] = (180./GTConst::Pi) * atan(pow(pow(dzdx->co[r][c], 2.0) + pow(dzdy->co[r][c], 2.0), 0.5));
				M[r][c] = (180./GTConst::Pi) * atan(pow(pow(dzdx[r][c], 2.0) + pow(dzdy[r][c], 2.0), 0.5));
			}else {
				M[r][c] = (double)undef;
			}
		}
	}

//	return M;
}


/////////////////////////////////////////////////////////////////////////////
// These routines coming from sparse_matrix_library.. S.C. 07.12.2013
/////////////////////////////////////////////////////////////////////////////


void product_matrix_using_lower_part_by_vector_plus_vector(double k, GeoVector<double>& out, const GeoVector<double>& y, const GeoVector<double>& x, const GeoVector<long>&  Li, const GeoVector<long>& Lp, GeoVector<double>& Lx){

	//calculates k*(y + Ax), where k is coefficient, y and x vectors, and A a SPD matrix defined with its lower diagonal part

	long i;

	product_using_only_strict_lower_diagonal_part(out, x, Li, Lp, Lx);

	for(i=1;i<x.size();i++){
#ifdef VERBOSE		
		printf("product_matrix_using_lower_part_by_vector_plus_vector-> i:%ld k:%e out:%e y:%e out:%e\n",i,k,out[i],y[i],k * (out[i] + y[i]));
#endif		
		out[i] = k * (out[i] + y[i]);
	}
}



//======

void product_using_only_strict_lower_diagonal_part(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<long>& Li, const GeoVector<long>& Lp, GeoVector<double>& Lx){
	long c, r ;

	for(size_t i=1;i<x.size();i++)
	{
	  product[i] = 0.0;
	}
	
	c = 1;
	for(size_t i=1;i<Li.size();i++){
		
		r = Li[i];

		if(r > c){


			product[c] += Lx[i] * (x[r] - x[c]);
			product[r] += Lx[i] * (x[c] - x[r]);
			
#ifdef VERBOSE				
		//printf("c:%ld -> %e i:%ld Lx:%e r:%ld %e c:%ld %e -> %e \n",c,product[c],i,Lx[i],r,x[r],c,x[c],Lx[i] * (x[r] - x[c]));
		//printf("r:%ld -> %e i:%ld Lx:%e c:%ld %e r:%ld %e -> %e \n",r,product[r],i,Lx[i],c,x[c],r,x[r],Lx[i] * (x[c] - x[r]));
#endif		   
			
		}else if(r < c){
#ifdef VERBOSE
			printf("product_using_only_strict_lower_diagonal_part r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Lp[c],Li.size(),Lp[x.size()-1],x.size());
#endif
            t_error("product_using_only_strict_lower_diagonal_part:matrix is not L, see function: " + std::string(__FUNCTION__));
		}


		if(i<Li.size()){
			while(i >= Lp[c]) c++;
		}
	}
}



//====================================

long BiCGSTAB_strict_lower_matrix_plus_identity_by_vector(double tol_rel, double tol_min, double tol_max, GeoVector<double>& x,
	const GeoVector<double>& b, GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx){


	//solve sistem (A+Iy)*x = B, find x
	//A M-matrix described by its lower diagonal part


	GeoVector<double> r0, r, p, v, s, t, diag, udiag, yy, z;
	double rho, rho1, alpha, omeg, beta, norm_r0;
	long i=0, j, maxiter;
	short sux;


	/*r0.resize(x.size() + 1);
	r.resize(x.size() + 1);
	p.resize(x.size() + 1);
	v.resize(x.size() + 1);
	s.resize(x.size() + 1);
	t.resize(x.size() + 1);
	diag.resize(x.size() + 1);
	udiag.resize(x.size() + 1);
	yy.resize(x.size() + 1);
	z.resize(x.size() + 1);*/

	r0.resize(x.size());
	r.resize(x.size());
	p.resize(x.size());
	v.resize(x.size());
	s.resize(x.size());
	t.resize(x.size());
	diag.resize(x.size());
	udiag.resize(x.size() -1 );
	yy.resize(x.size());
	z.resize(x.size());
	
	
	
	get_diag_strict_lower_matrix_plus_identity_by_vector(diag, udiag, y, Li, Lp, Lx);

	//product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(r, x, y, Li, Lp, Lx);

	///printf("==== >> x->nl=%ld x->nh=%ld ",x->nl,x->nh);
	//stop_execution();

	
	for (j=GTConst::nl;j<x.size();j++ ) {
   		
    	r[j] = b[j];
		r0[j] = r[j];
		p[j] = 0.;
		v[j] = 0.;
    }


	norm_r0 = norm_2(r0, GTConst::nl, r0.size());
#ifdef VERBOSE
	printf("BiCGSTAB_strict norm_r0: %f\n",norm_r0);
#endif 

	rho = 1.;
	alpha = 1.;
	omeg = 1.;

	maxiter = (long)(x.size()/100.); 
	if (maxiter < 100) maxiter = 100;
	
  
	while ( i<maxiter && norm_2(r, GTConst::nl, r.size()) > Fmax( tol_min , Fmin( tol_max , tol_rel*norm_r0) ) ) {
		
		rho1 = product(r0, r);

		beta = (rho1/rho)*(alpha/omeg);

		rho = rho1;

	
		for (j=GTConst::nl;j<x.size();j++ ) {
			
			p[j] = r[j] + beta*(p[j] - omeg*v[j]);
		}

		sux=tridiag(0, 0, 0, x.size(), udiag, diag, udiag, p, yy);
		
		if(sux==0) return(-1);

		product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(v, yy, y, Li, Lp, Lx);

		alpha = rho/product(r0, v);

		for (j=GTConst::nl;j<x.size();j++ ) {
			s[j] = r[j] - alpha*v[j];
		}

		if(norm_2(s,GTConst::nl,s.size())>1.E-10){
			sux=tridiag(0, 0, 0, x.size(), udiag, diag, udiag, s, z);
			if(sux==0) return(-1);

			product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(t, z, y, Li, Lp, Lx);

			omeg = product(t, s)/product(t, t);

			for (j=GTConst::nl;j<x.size();j++ ) {

				x[j] += (alpha*yy[j] + omeg*z[j]);
				r[j] = s[j] - omeg*t[j];
			}

		}else{


			for (j=GTConst::nl;j<x.size();j++ ) {
				x[j] += alpha*yy[j];
				r[j] = s[j];
			}
		}
		
		i++;
	}

	//free_doublevector(r0);
	//free_doublevector(r);
	//free_doublevector(p);
	//free_doublevector(v);
	//free_doublevector(s);
	//free_doublevector(t);
	//free_doublevector(diag);
	//free_doublevector(udiag);
	//free_doublevector(yy);
	//free_doublevector(z);

	return i;

}



//===============================

short tridiag(short a, long r, long c, long nx, const GeoVector<double>& diag_inf, const GeoVector<double>& diag, const GeoVector<double>& diag_sup, const GeoVector<double>& b, GeoVector<double>& e)
{

long j;
double bet;

GeoVector<double> gam;


	gam.resize(nx);

	
if(diag[1]==0.0){
	printf("type=%d r=%ld c=%ld\n",a,r,c);
	t_error("Error 1 in tridiag");
}


	bet=diag[1];
	e[1]=b[1]/bet;

for(j=2;j<nx;j++){
	gam[j]=diag_sup[j-1]/bet;
	bet=diag[j]-diag_inf[j-1]*gam[j];
	if(bet==0.0){
		printf("type=%d r=%ld c=%ld\n",a,r,c);
		printf("l=%ld diag(l)=%f diag_inf(l-1)=%f diag_sup(l-1)=%f\n",j,diag[j],diag_inf[j-1],diag_sup[j-1]);
		printf("Error 2 in tridiag\n");
		return 0;
	}

	e[j]=(b[j]-diag_inf[j-1]*e[j-1])/bet;
}

//Backsubstitution
for(j=(nx-2);j>=1;j--){
	e[j]-=gam[j+1]*e[j+1];
}

//free_doublevector(gam);

return 1;

}

//void get_diag_strict_lower_matrix_plus_identity_by_vector(DOUBLEVECTOR *diag, DOUBLEVECTOR *udiag, DOUBLEVECTOR *y,
//	LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){



void get_diag_strict_lower_matrix_plus_identity_by_vector(GeoVector<double>& diag,  GeoVector<double>& udiag, const GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx){


	long i, r, c;
	
	//find diagonal and upper diagonal of matrix A+Iy, where A is described by its strict lower diagonal part

	for(i=1;i<diag.size();i++){
		diag[i] = y[i];
		if(i<diag.size()-1) udiag[i] = 0.;
	}

	c = 1;
	for(i=1;i<Li.size();i++){
		r = Li[i];

		diag[c] -= Lx[i];
		diag[r] -= Lx[i];
		
		if(r == c+1) udiag[c] = Lx[i];
		if(i<Li.size()-1){
			while(i >= Lp[c]) c++;
		}
	}
}



//void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(DOUBLEVECTOR *product, DOUBLEVECTOR *x, DOUBLEVECTOR *y,
//																	LONGVECTOR *Li, LONGVECTOR *Lp, DOUBLEVECTOR *Lx){


void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<double>& y,
																	const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx){

	long i, r, c;

	//calculate (A+Iy)*x, A described by its strict lower diagonal part

//	for(i=1;i<=x->nh;i++){
	for(i=1;i<x.size();i++){
	//  product->co[i] = y->co[i] * x->co[i];
		product[i] = y[i] * x[i];
	}

	c = 1;

	for(i=1;i<Li.size();i++){
		r = Li[i];

		if(r > c){
			//printf("product_using_only_strict_lower_diagonal_part: r:%ld c:%ld i:%ld Lx:%f x->co[c]:%f xi->co[r]: %f\n",r,c,i,Lx[i],x[c],x[r]);

						product[c] += Lx[i] * (x[r] - x[c]);
						product[r] += Lx[i] * (x[c] - x[r]);
		}else if(r < c){
			//printf(" product_using_only_strict_lower_diagonal_part_plus_identity_by_vector r:%ld c:%ld i:%ld Lp[c]:%ld tot:%ld %ld %ld\n",r,c,i,Lp[c],Li.size(),Lp[x.size()],x.size());
            t_error("product_using_only_strict_lower_diagonal_part_plus_identity_by_vector: matrix is not L, see function: " + std::string(__FUNCTION__));
		}

		if(i<Li.size()){
			while(i >= Lp[c]) c++;
		}
	}
}

double product(const GeoVector<double>& a, const GeoVector<double>& b){

	double p=0.;

	long i ;

	for(i=1;i<a.size();i++){
		p += a[i] * b[i];
	}

	return(p);
}
