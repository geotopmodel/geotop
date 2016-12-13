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
    if(p != parent && mkdirp(parent, mode) != 0) {
#if defined(__CYGWIN__)
        return 0;
#else
        return -1;
#endif
    }

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
    
	int ret = mkdirp(basedir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1 == ret) {
		fprintf(stderr, "ERROR: Unable to create parent directory `%s` for name `%s`. Exiting.\n", basedir, name);
		exit(1);
	}
	free(lBasePtr); //lBasePtr cannot be freed as long basedir might be in use - see 'man 3 dirname'

	if ((fp=fopen(name,mode))==NULL) {
		printf("%s", name);    
		t_error(" The specified file could not be opened ");
		
		return NULL;
	} else {
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
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


//----------------------------
void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n){

    long i ;
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
}

//----------------------------------------

void find_min_max(GeoMatrix<double>& M, long novalue, double *max, double *min){

	long r, c, nr=M.getRows(), nc=M.getCols();

	*max=-1.E99;
	*min=1.E99;

	for (r=1; r<nr; r++) {
		for (c=1; c<nc; c++) {
			if ((long)M[r][c] != novalue) {
				if(*max < M[r][c]) *max = M[r][c];
				if(*min > M[r][c]) *min = M[r][c];
			}
		}
	}
}

//-----------------------------------------
double norm_1(const GeoVector<double>& V, long nbeg, long nend){

	long l;
	double N=0.0;
	for(l=nbeg;l<=nend;l++){
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

// adapted to GEOMatrix/Vector to get rid of fluidturtle datastructure: to double check if things are ok SC 13.12.2016

void sky_view_factor(GeoMatrix<double>& sky, long N, TInit *UV, GeoMatrix<double>& input, GeoMatrix<short>& convess, long novalue)
{
    long i,j,t,m,n,p,q,h,k,r,s; //counters
    double deltateta; //amplitude of the angles in which the horizon is divided
    GeoMatrix<double> alfa; //matrices with the angles of the direction
    GeoVector<double> v; //vector with the view factor of the current pixel for one of the N parts
    GeoVector<double> vv; //vector with the minimum view factor of the current pixel for one of the N parts
    double vvv; //mean of the sky view for a pixel of the N parts
    long nr=input.getRows()-1;
    long nc=input.getCols()-1;
    
    if(sky.getRows()-1!=nr) t_error("Sky view factor fatal error, number of rows not consistent");
    if(sky.getCols()-1!=nc) t_error("Sky view factor fatal error, number of cols not consistent");
    
    // Computation of the matrix with the angles of the direction
    // alfa=new_doublematrix(2*nr-1,2*nc-1);
    // ALLOCATION: please note that I allocate 1 element more for columns and row.. " TO  double checks 
    alfa.resize(2*nr,2*nc,-9999.0);
    

    for(i=1;i<=2*nr-1;i++){
        //  for(j=1;j<=2*input->nch-1;j++){
        for(j=1;j<=2*nc-1;j++){
            //  if(i<=input->nrh && j<input->nch){
            if(i<=nr && j<nc){
                //  alfa->co[i][j]=3.0/2.0*GTConst::Pi+atan(((input->nrh-i)*UV->U->co[1])/((input->nch-j)*UV->U->co[1]));
                alfa[i][j]=3.0/2.0*GTConst::Pi+atan(((nr-i)*UV->U[1])/((nc-j)*UV->U[1]));
            }
            // if(i>input->nrh && j<=input->nch){
            if(i>nr && j<=nc){
                //  alfa->co[i][j]=GTConst::Pi+atan(((input->nch-j)*UV->U->co[1])/((i-input->nrh)*UV->U->co[1]));
                alfa[i][j]=GTConst::Pi+atan(((nc-j)*UV->U[1])/((i-nr)*UV->U[1]));
            }
            // if(i>=input->nrh && j>input->nch){
            if(i>=nr && j>nc){
                //  alfa->co[i][j]=GTConst::Pi/2.0+atan(((i-input->nrh)*UV->U->co[1])/
                //                        ((j-input->nch)*UV->U->co[1]));
                alfa[i][j]=GTConst::Pi/2.0+atan(((i-nr)*UV->U[1])/
                                                ((j-nc)*UV->U[1]));
                
            }
            //  if(i<input->nrh && j>=input->nch){
            if(i<nr && j>=nc){
                //   alfa->co[i][j]=atan(((j-input->nch)*UV->U->co[1])/((input->nrh-i)*UV->U->co[1]));
                alfa[i][j]=atan(((j-nc)*UV->U[1])/((nr-i)*UV->U[1]));
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
    
    // to carefully checking this below..
    
    v.resize(N+1);
    vv.resize(N+1);
    deltateta=2.0*GTConst::Pi/N;
    
    // for(i=1;i<=input->nrh;i++){
    for(i=1;i<=nr;i++){
        //	for(j=1;j<=input->nch;j++){
        for(j=1;j<=nc;j++){
            //	if ((long)input->co[i][j]!=novalue){ //computation only of novalue pixels
            if ((long)input[i][j]!=novalue){ //computation only of novalue pixels
                for(t=1;t<=N;t++){
                    v[t]=1.0;
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
                            if (alfa[h][k]>=(t-1)*deltateta && alfa[h][k]<t*deltateta){
                                r=h-m+1;
                                s=k-n+1;
                                // if (convess->co[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
                                if (convess[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
                                    //		                        vv->co[t]=1-sin(atan((input->co[r][s]-input->co[i][j])
                                    //		                         	          /(sqrt(pow((r-i),2)+pow((s-j),2))*UV->U->co[1])));
                                    vv[t]=1-sin(atan((input[r][s]-input[i][j])
                                                     /(sqrt(pow((r-i),2)+pow((s-j),2))*UV->U[1])));
                                    if (vv[t]<v[t]){
                                        v[t]=vv[t];
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                vvv=0.0;
                for(t=1;t<=N;t++){
                    vvv=vvv+v[t];
                }
                //  sky->co[i][j]=(1.0/N*vvv);
                sky[i][j]=(1.0/N*vvv);
            }
        }
        // 	printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",100.0*(double)i/(double)sky->nrh);
        printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",100.0*(double)i/(double)sky.getRows()-1);
    }
    
    //#free_doublematrix(alfa);
    //free_doublevector(v);
    //free_doublevector(vv);
    
}

//***************************************************************************

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
// overloade function noori
void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n){

	long r, c, nr, nc, ir, ic, i;
//	DOUBLEVECTOR *values;
    GeoVector<double> values;
    long cnt;

	values.resize((2*n+1)*(2*n+1));

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
								values[cnt]=Zin[r+ir][c+ic];
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
					Zout[r][c] += values[i]/(double)cnt;
				}

			}else {

			//	Zout->co[r][c] = (double)novalue;
				Zout[r][c] = (double)novalue;

			}

		}
	}

//	free_doublevector(values);
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



//overloaded function noori
void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M){

	long r, c;
	long nc=topo.getCols();
	long nr=topo.getRows();

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

	size_t i;

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
			
#ifdef VERY_VERBOSE				
		printf("c:%ld -> %e i:%ld Lx:%e r:%ld %e c:%ld %e -> %e \n",c,product[c],i,Lx[i],r,x[r],c,x[c],Lx[i] * (x[r] - x[c]));
		printf("r:%ld -> %e i:%ld Lx:%e c:%ld %e r:%ld %e -> %e \n",r,product[r],i,Lx[i],c,x[c],r,x[r],Lx[i] * (x[c] - x[r]));
#endif		   
			
		}else if(r < c){
#ifdef VERBOSE
			printf("product_using_only_strict_lower_diagonal_part r:%ld c:%ld i:%ld Ap[c]:%ld tot:%ld %ld %ld\n",r,c,i,Lp[c],Li.size(),Lp[x.size()-1],x.size());
#endif
            t_error("product_using_only_strict_lower_diagonal_part:matrix is not L, see function: " + std::string(__FUNCTION__));
		}


 		if(i < Li.size() && (long)Li.size() > 0){
			while(c < (long)Lp.size() && (long)i >= Lp[c])
                c++;
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
	long i=0, maxiter;
	short sux;
	size_t j;


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
    if (norm_r0 != norm_r0) {
#ifdef VERBOSE
        printf(" BiCGSTAB_strict norm_r0: %f\n",norm_r0);
#endif
        t_error("Fatal Error! NAN in a variable. See failing report.");

        
    }
#ifdef VERBOSE
	printf(" BiCGSTAB_strict norm_r0: %f\n",norm_r0);
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


void get_diag_strict_lower_matrix_plus_identity_by_vector(GeoVector<double>& diag,  GeoVector<double>& udiag, const GeoVector<double>& y, const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx){


	long r, c;
	size_t i;
	
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
		if(i<Li.size()-1 && (long)i > 0){
			while(c < (long)Lp.size() && (long)i >= Lp[c])
                c++;
		}
	}
}


void product_using_only_strict_lower_diagonal_part_plus_identity_by_vector(GeoVector<double>& product, const GeoVector<double>& x, const GeoVector<double>& y,
																	const GeoVector<long>& Li, const GeoVector<long>& Lp, const GeoVector<double>& Lx){

	long r, c;
  	size_t i;

	//calculate (A+Iy)*x, A described by its strict lower diagonal part

	for(i=1;i<x.size();i++){
		product[i] = y[i] * x[i];
	}

	c = 1;

	for(i=1;i<Li.size();i++){
		r = Li[i];

		if(r > c){
						product[c] += Lx[i] * (x[r] - x[c]);
						product[r] += Lx[i] * (x[c] - x[r]);
		}else if(r < c){
            t_error("product_using_only_strict_lower_diagonal_part_plus_identity_by_vector: matrix is not L, see function: " + std::string(__FUNCTION__));
		}

		if(i<Li.size() && (long)i > 0){
			while(c < (long)Lp.size() && (long)i >= Lp[c])
                c++;
		}
	}
}

double product(const GeoVector<double>& a, const GeoVector<double>& b){

	double p=0.;

	size_t i;

	for(i=1;i<a.size();i++){
		p += a[i] * b[i];
	}

	return(p);
}
