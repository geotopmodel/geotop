
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi 
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#include "tables.h"
		
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_activelayerdepth_up(long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long nmax=sl->pa->nch;
	long l;//counter
	short out=0;
	
	n = nmax;
	
	if(sl->SS->T->co[n][i]>=thresh){
		for(l=1;l<=n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
	
	}else{
	
		do{
			n--;
			if(n==1) out=-1;
			if(sl->SS->T->co[n+1][i]<thresh && sl->SS->T->co[n][i]>=thresh) out=1;
		}while(out==0);
		
		if(out==1){
			
			for(l=1;l<=n;l++){
				table += sl->pa->co[ty][jdz][l];
			}
			
			table += (1.-sl->SS->thi->co[n+1][i]/(sl->SS->thi->co[n+1][i]+sl->th->co[n+1][i]-sl->pa->co[ty][jres][n+1]))*sl->pa->co[ty][jdz][n+1];
		}
	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_activelayerdepth_dw(long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long nmax=sl->pa->nch;
	long l;//counter
	short out=0;
	
	n = 1;
	
	if(sl->SS->T->co[n][i]<thresh){
		
		do{
			n++;
			if(n==nmax) out=-1;
			if(sl->SS->T->co[n][i]>=thresh && sl->SS->T->co[n-1][i]<thresh) out=1;
		}while(out==0);
		
		for(l=1;l<n-1;l++){
			table += sl->pa->co[ty][jdz][l];
		}
		
		if(out==1){
			table += (sl->SS->thi->co[n-1][i]/(sl->SS->thi->co[n-1][i]+sl->th->co[n-1][i]-sl->pa->co[ty][jres][n-1]))*sl->pa->co[ty][jdz][n-1];
		}else {
			table = 0.;
		}

	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_watertabledepth_up(long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long nmax=sl->pa->nch;
	long l;//counter
	short out=0;
	
	n = nmax;
	
	if(sl->Ptot->co[n][i]<thresh){
		for(l=1;l<=n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
		
	}else{
		
		do{
			n--;
			if(n==1) out=-1;
			if(sl->Ptot->co[n+1][i]>=thresh && sl->Ptot->co[n][i]<thresh) out=1;
		}while(out==0);
		
		if(out==1){
			
			for(l=1;l<=n;l++){
				table += sl->pa->co[ty][jdz][l];
			}
			
			table += ( 0.5*sl->pa->co[ty][jdz][n+1] - sl->Ptot->co[n+1][i]*0.5*(sl->pa->co[ty][jdz][n]+sl->pa->co[ty][jdz][n+1])/(sl->Ptot->co[n+1][i]-sl->Ptot->co[n][i]) );
		}
	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_watertabledepth_dw(long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long nmax=sl->pa->nch;
	long l;//counter
	short out=0;
	
	n = 1;
		
	if(sl->Ptot->co[n][i]<thresh){
		
		do{
			n++;
			if(n==nmax) out=-1;
			if(sl->Ptot->co[n][i]>=thresh && sl->Ptot->co[n-1][i]<thresh) out=1;
		}while(out==0);
							
		for(l=1;l<n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
		
		if (out==1){
			table += ( 0.5*sl->pa->co[ty][jdz][n] - sl->Ptot->co[n][i]*0.5*(sl->pa->co[ty][jdz][n-1]+sl->pa->co[ty][jdz][n])/(sl->Ptot->co[n][i]-sl->Ptot->co[n-1][i]) );
		}else {
			table += sl->pa->co[ty][jdz][n];
		}
	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
