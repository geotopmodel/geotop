
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

#include "constant.h"
#include "keywords_file.h"
#include "struct.geotop.h"
#include "tables.h"

extern T_INIT *UV;
		
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_activelayerdepth(long r, long c, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long l;//counter
	long ty=sl->type->co[r][c];
	short out=0;
	
	n = sl->pa->nch;
	
	if(sl->T->co[n][r][c]>=thresh){
		for(l=1;l<=n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
	
	}else{
	
		do{
			n--;
			if(n==1) out=-1;
			if(sl->T->co[n+1][r][c]<thresh && sl->T->co[n][r][c]>=thresh) out=1;
		}while(out==0);
		
		if(out==1){
			
			for(l=1;l<=n;l++){
				table += sl->pa->co[ty][jdz][l];
			}
			
			table += (1.-sl->thice->co[n+1][r][c]/(sl->thice->co[n+1][r][c]+sl->th->co[n+1][r][c]-sl->pa->co[ty][jres][n+1]))*sl->pa->co[ty][jdz][n+1];
		}
	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double find_watertabledepth(long r, long c, SOIL *sl)

{
	double table=0.0;	
	double thresh=0.0;
	long n;// number of layer below the threshold
	long l;//counter
	long ty=sl->type->co[r][c];
	short out=0;
	
	n = sl->pa->nch;
	
	if(sl->Ptot->co[n][r][c]<thresh){
		for(l=1;l<=n;l++){
			table += sl->pa->co[ty][jdz][l];
		}
		
	}else{
		
		do{
			n--;
			if(n==1) out=-1;
			if(sl->Ptot->co[n+1][r][c]>=thresh && sl->Ptot->co[n][r][c]<thresh) out=1;
		}while(out==0);
		
		if(out==1){
			
			for(l=1;l<=n;l++){
				table += sl->pa->co[ty][jdz][l];
			}
			
			table += ( 0.5*sl->pa->co[ty][jdz][n+1] - sl->Ptot->co[n+1][r][c]*0.5*(sl->pa->co[ty][jdz][n]+sl->pa->co[ty][jdz][n+1])
					  /(sl->Ptot->co[n+1][r][c]-sl->Ptot->co[n][r][c]) );
		}
	}
	
	return table;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
