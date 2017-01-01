/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */


#include <string>
#include "geomorphology.h"
#include "read_command_line.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <iostream>
#include <meteoio/MeteoIO.h>


/*WORKING_POSITION=SEEK_SET;*/

char  EXTERNAL_FILE_NAME[256],OLD_NAME[256];

long int EXTERNAL_FILE_POSITION=SEEK_SET;


short OPENYES=0;

long IO_FILES_COUNTER=0;
long IO_STRINGS_COUNTER=0;
long  WORKING_POSITION=0;
long IO_PARMS_COUNTER=0;


//----------------------------
void multipass_topofilter(long ntimes, GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n){
    
    long i ;
    GeoMatrix<double> M;
    
    Zout=Zin;
    
    for (i=1; i<=ntimes; i++) {
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
    // ALLOCATION: please note that I allocate 1 element more for columns and row.. " TODO  double checks
    alfa.resize(2*nr,2*nc,-9999.0);
    
    
    for(i=1;i<=2*nr-1;i++){
        for(j=1;j<=2*nc-1;j++){
              if(i<=nr && j<nc){
                alfa[i][j]=3.0/2.0*GTConst::Pi+atan(((nr-i)*UV->U[1])/((nc-j)*UV->U[1]));
            }
            if(i>nr && j<=nc){
            
                alfa[i][j]=GTConst::Pi+atan(((nc-j)*UV->U[1])/((i-nr)*UV->U[1]));
            }
                 if(i>=nr && j>nc){
                alfa[i][j]=GTConst::Pi/2.0+atan(((i-nr)*UV->U[1])/
                                                ((j-nc)*UV->U[1]));
                
            }
            if(i<nr && j>=nc){
                 alfa[i][j]=atan(((j-nc)*UV->U[1])/((nr-i)*UV->U[1]));
            }
        }
    }
    
    // Computation of matrix with sky view factor:
     for(i=1;i< sky.getRows();i++){

        for(j=1;j<sky.getCols();j++){
            sky[i][j]=(double)novalue;
        }
    }
    
    // to carefully checking this below..
    
    v.resize(N+1);
    vv.resize(N+1);
    deltateta=2.0*GTConst::Pi/N;
    
    for(i=1;i<=nr;i++){
 
        for(j=1;j<=nc;j++){
 
            if ((long)input[i][j]!=novalue){ //computation only of novalue pixels
                for(t=1;t<=N;t++){
                    v[t]=1.0;
                }
                m=nr-i+1;
                n=nc-j+1;
                p=m+nr-1;
                q=n+nc-1;
                for(h=m;h<=p;h++){
                    for(k=n;k<=q;k++){
                        for(t=1;t<=N;t++){
                            if (alfa[h][k]>=(t-1)*deltateta && alfa[h][k]<t*deltateta){
                                r=h-m+1;
                                s=k-n+1;
                                if (convess[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
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
                sky[i][j]=(1.0/double(N)*vvv);
            }
        }
         printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",100.0*(double)i/(double)sky.getRows()-1);
    }
    
    
}

//***************************************************************************

void topofilter(GeoMatrix<double>& Zin, GeoMatrix<double>& Zout, long novalue, long n)
{
    
    long r, c, nr, nc, ir, ic, i;
    GeoVector<double> values;
    long cnt;
    
    values.resize((2*n+1)*(2*n+1)+1);
    // the +1 above is needed to prevent a memory corruptions on glibc on Linux boxes (S.C. 16.12.2016)  
    nr=Zin.getRows();
    nc=Zin.getCols();
    
    for (r=1; r<nr; r++) {
        for (c=1; c<nc; c++) {
            if ((long)Zin[r][c] != novalue) {
                cnt=0;
                for (ir=-n; ir<=n; ir++) {
                    for (ic=-n; ic<=n; ic++) {
                        if (r+ir>=1 && r+ir<=nr && c+ic>=1 && c+ic<=nc) {
                            if((long)Zin[r+ir][c+ic]!=novalue){
                                cnt++;
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
                
                Zout[r][c] = 0.;
                for (i=1; i<=cnt; i++) {
                    Zout[r][c] += values[i]/(double)cnt;
                }
                
            }else {
                
                Zout[r][c] = (double)novalue;
                
            }
            
        }
    }
    
}

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

void find_max_slope(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M)
{
    
    long r, c;
    long nc=topo.getCols();
    long nr=topo.getRows();
    
    
    M.resize(nr, nc);
    for(r=1;r<nr;r++){
        for(c=1;c<nc;c++){
            if((long)topo[r][c]!=undef){
                M[r][c] = (180./GTConst::Pi) * atan(pow(pow(dzdx[r][c], 2.0) + pow(dzdy[r][c], 2.0), 0.5));
            }else {
                M[r][c] = (double)undef;
            }
        }
    }
}

//overloaded function noori
void find_aspect(GeoMatrix<double>& topo, GeoMatrix<double>& dzdx, GeoMatrix<double>& dzdy, long undef, GeoMatrix<double>& M){

	long r, c;
	long nc=topo.getCols();
	long nr=topo.getRows();
	M.resize(nr,nc,-9999.0);
	for(r=1;r<=nr;r++){
	for(r=1;r<nr;r++){
		for(c=1;c<nc;c++){
			if((long)topo[r][c]!=undef){
				M[r][c] = (180./GTConst::Pi) * (3.0 / 2.0 * GTConst::Pi - atan2(dzdy[r][c], dzdx[r][c]));
				if (M[r][c]>=360.0) M[r][c] -= 360.0;
			}else {
				M[r][c] = (double)undef;
			}
		}
	}
}




}
