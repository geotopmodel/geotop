#include "turtle.h"
#include "dtm_resolution.h"

//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------

//1 known, 2 unknown
void reduce_resolution(long n, DOUBLEMATRIX *DTM1, DOUBLEMATRIX *DTM2, T_INIT *UV1, T_INIT *UV2){

	long r,c,rtot,ctot,rr,cc,cont,Rr=0,Rc=0;
	long *contvt;
	
	rtot=(long)(DTM1->nrh/(double)n);
	ctot=(long)(DTM1->nch/(double)n);
	for(r=1;r<=rtot;r++){
		for(c=1;c<=ctot;c++){
			DTM2->co[r][c]=0.0;
			cont=0;
			for(rr=1;rr<=n;rr++){
				for(cc=1;cc<=n;cc++){
					if(DTM1->co[n*(r-1)+rr][n*(c-1)+cc]!=UV1->V->co[2]){
						DTM2->co[r][c]+=DTM1->co[n*(r-1)+rr][n*(c-1)+cc];
						cont+=1;
					}
				}
			}
			if(cont==0){
				DTM2->co[r][c]=UV1->V->co[2];
			}else{
				DTM2->co[r][c]/=(double)cont;
			}
		}
	}
	
	if(n*rtot<DTM1->nrh){
		contvt=(long *)malloc(ctot*sizeof(long));
		cont=0;
		for(c=1;c<=ctot;c++){
			DTM2->co[rtot+1][c]=0.0;
			contvt[c-1]=0;	
		}
		do{
			if(n*rtot+cont+1<=DTM1->nrh){
				cont+=1;
				for(c=1;c<=ctot;c++){
					for(cc=1;cc<=n;cc++){
						if(DTM1->co[n*rtot+cont][n*(c-1)+cc]!=UV1->V->co[2]){
							DTM2->co[rtot+1][c]+=DTM1->co[n*rtot+cont][n*(c-1)+cc];
							contvt[c-1]+=1;
						}
					}
				}
			}
		}while(n*rtot+cont<DTM1->nrh);
		Rr=cont;
		for(c=1;c<=ctot;c++){
			if(contvt[c-1]==0){
				DTM2->co[rtot+1][c]=UV1->V->co[2];
			}else{
				DTM2->co[rtot+1][c]/=(double)contvt[c-1];	
			}
		}
		free(contvt);
	}

	if(n*ctot<DTM1->nch){
		contvt=(long *)malloc(rtot*sizeof(long));
		cont=0;
		for(r=1;r<=rtot;r++){
			DTM2->co[r][ctot+1]=0.0;
			contvt[r-1]=0;
		}
		do{
			if(n*ctot+cont+1<=DTM1->nch){
				cont+=1;				
				for(r=1;r<=rtot;r++){
					for(rr=1;rr<=n;rr++){
						if(DTM1->co[n*(r-1)+rr][n*ctot+cont]!=UV1->V->co[2]){
							DTM2->co[r][ctot+1]+=DTM1->co[n*(r-1)+rr][n*ctot+cont];
							contvt[r-1]+=1;
						}
					}
				}
			}
		}while(n*ctot+cont<DTM1->nch);
		Rc=cont;
		for(r=1;r<=rtot;r++){
			if(contvt[r-1]==0){
				DTM2->co[r][ctot+1]=UV1->V->co[2];
			}else{
				DTM2->co[r][ctot+1]/=(double)contvt[r-1];
			}
		}
		free(contvt);
	}
	
	if(Rr!=0 && Rc!=0){
		DTM2->co[rtot+1][ctot+1]=0.0;
		cont=0;
		for(rr=n*rtot+1;rr<=n*rtot+Rr;rr++){
			for(cc=n*ctot+1;cc<=n*ctot+Rc;cc++){
				if(DTM1->co[rr][cc]!=UV1->V->co[2]){
					DTM2->co[rtot+1][ctot+1]+=DTM1->co[rr][cc];
					cont+=1;
				}
			}
		}
		if(cont==0){
			DTM2->co[rtot+1][ctot+1]=UV1->V->co[2];
		}else{
			DTM2->co[rtot+1][ctot+1]/=(double)cont;
		}
	}

	UV2->U->co[1]=(UV1->U->co[1]*DTM1->nrh)/(double)DTM2->nrh;
	UV2->U->co[2]=(UV1->U->co[2]*DTM1->nch)/(double)DTM2->nch;
	UV2->U->co[3]=UV1->U->co[3];
	UV2->U->co[4]=UV1->U->co[4];
	UV2->V->co[1]=UV1->V->co[1];	
	UV2->V->co[2]=UV1->V->co[2];	
	
}
		
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------

//1 known, 2 unknown
void amplify_resolution(long n, long Rr, long Rc, DOUBLEMATRIX *DTM1, DOUBLEMATRIX *DTM2, T_INIT *UV1, T_INIT *UV2){

	long r,c,rr,cc,rtot,ctot;
	
	rtot=DTM1->nrh;
	ctot=DTM1->nch;
	if(Rr>0) rtot-=1;
	if(Rc>0) ctot-=1;
	
	for(r=1;r<=rtot;r++){
		for(c=1;c<=ctot;c++){
			for(rr=1;rr<=n;rr++){
				for(cc=1;cc<=n;cc++){
					DTM2->co[rr+(r-1)*n][cc+(c-1)*n]=DTM1->co[r][c];
				}
			}
		}
	}
	
	if(Rr>0){
		for(c=1;c<=ctot;c++){
			for(rr=1;rr<=Rr;rr++){
				for(cc=1;cc<=n;cc++){
					DTM2->co[n*rtot+rr][cc+(c-1)*n]=DTM1->co[rtot+1][c];
				}
			}
		}
	}
	
	if(Rc>0){
		for(r=1;r<=rtot;r++){
			for(rr=1;rr<=n;rr++){
				for(cc=1;cc<=Rc;cc++){
					DTM2->co[rr+(r-1)*n][n*ctot+cc]=DTM1->co[r][ctot+1];
				}
			}
		}
	}
		
	if(Rr>0 && Rc>0){
		for(rr=1;rr<=Rr;rr++){
			for(cc=1;cc<=Rc;cc++){
				DTM2->co[n*rtot+rr][n*ctot+cc]=DTM1->co[rtot+1][ctot+1];
			}
		}
	}
	
	UV2->U->co[1]=(UV1->U->co[1]*DTM1->nrh)/(double)DTM2->nrh;
	UV2->U->co[2]=(UV1->U->co[2]*DTM1->nch)/(double)DTM2->nch;
	UV2->U->co[3]=UV1->U->co[3];
	UV2->U->co[4]=UV1->U->co[4];
	UV2->V->co[1]=UV1->V->co[1];	
	UV2->V->co[2]=UV1->V->co[2];	
	
}

//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------

