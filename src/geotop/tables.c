
/* STATEMENT:
 
 Geotop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 Geotop 2.0.0 - 31 Oct 2013
 
 Copyright (c), 2013 - Stefano Endrizzi 
 
 This file is part of Geotop 2.0.0
 
 Geotop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 Geotop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the Geotop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */

#include "constants.h"
#include "struct.geotop.h"
#include "tables.h"

extern T_INIT *UV;
		
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
		
	if(sl->SS->T->co[n][i]<thresh && nmax>1){
		
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

double find_watertabledepth_up(double Z, long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long l;//counter
	short out=0;
		
	n = nlayer(Z, sl->pa->co[ty][jdz], sl->pa->nch, -1);

	if(n>1){
	
		if(sl->Ptot->co[n][i] < thresh){
			
			for(l=1;l<=n;l++){
				table += sl->pa->co[ty][jdz][l];
			}
			
		}else{
			
			do{
				n--;
				if(n==1) out=-1;
				if(sl->Ptot->co[n+1][i] >= thresh && sl->Ptot->co[n][i] < thresh) out=1;
			}while(out==0);
			
			if(out==1){
				
				for(l=1;l<=n;l++){
					table += sl->pa->co[ty][jdz][l];
				}
				
				table += 0.5 * sl->pa->co[ty][jdz][n+1];
				table -= 0.5 * (sl->pa->co[ty][jdz][n]+sl->pa->co[ty][jdz][n+1]) * (sl->Ptot->co[n+1][i]-thresh) / (sl->Ptot->co[n+1][i]-sl->Ptot->co[n][i]);
				
			}
		}
	}
	
	if(table>Z) table = Z;

	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_watertabledepth_dw(double Z, long i, long ty, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long nmax=sl->pa->nch;
	long l;//counter
	short out=0;
	
	n = 3;
		
	if(sl->Ptot->co[n][i] < thresh){
				
		do{
			n++;
			if (n>=nmax) out=-1;
			if (n<=nmax && sl->Ptot->co[n][i] >= thresh && sl->Ptot->co[n-1][i] < thresh) out=1;
			// Previously written:
//			if (n==nmax) out=-1;
//			if (sl->Ptot->co[n][i] >= thresh && sl->Ptot->co[n-1][i] < thresh) out=1;
		} while(out==0);
							
		for(l=1;l<n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
		
		if (out==1){
			table += ( 0.5*sl->pa->co[ty][jdz][n] - (sl->Ptot->co[n][i]-thresh)*0.5*(sl->pa->co[ty][jdz][n-1]+sl->pa->co[ty][jdz][n])/(sl->Ptot->co[n][i]-sl->Ptot->co[n-1][i]) );
		}else {
			table += sl->pa->co[ty][jdz][n];
		}
	}
	
	if(table>Z) table = Z;
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

//if d=1 finds the layer up, if d=-1 the layer down
long nlayer(double D, double *dz, long max, short d){
	
	double z=0.;
	long l=1;	
		
	do{ 
		if(l>1) z += dz[l-1]/2.;
		z += dz[l]/2.;
		l++;
	}while(z<D && l<=max);
	
		
	if (z>=D){
		if (d==-1){
			return(l-1);
		}else {
			if (fabs(z-D)<1.E-5 || l-1==1){
				return(l-1);
			}else {
				return(l-2);
			}
		}
	}else {
		return(l-1);
	}
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------


		
		
