#include  "turtle.h"
#include  "tensor3D.h"
/* Note that depth is the first indices and that the indices were pernutated
with respect to NR */

/*-----------------------------------------------------------------------*/


double    ***d3tensor( long nrl, long nrh, long ncl, long nch, long ndl, long ndh)

{

long i,j,nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
double ***t;


t=(double ***)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));

if(!t) t_error("Allocation failure of a double tensor pointers");
t+=NR_END;
t-=nrl;

t[nrl]=(double  **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
if(!t[nrl])  t_error("Allocation failure in a double tensors rows pointer" );
t[nrl]+=NR_END;
t[nrl]-=ncl;

t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
if(!t[nrl][ncl]) t_error("Allocation failure in a double tensor colunmns pointer");

t[nrl][ncl]+=NR_END;
t[nrl][ncl]-=ndl;


for(j=ncl+1;j<=nch; j++) t[nrl][j]=t[nrl][j-1]+ndep;
for(i=nrl+1;i<=nrh;i++){
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;


}

return t;


}


/*-----------------------------------------------------------------------*/


DOUBLETENSOR *new_doubletensor(long ndh,long nrh,long nch)


{

DOUBLETENSOR *m;

		  m=(DOUBLETENSOR *)malloc(sizeof(DOUBLETENSOR));
		  if (!m) t_error("allocation failure in new_doubletensor()");
		  m->isdynamic=isDynamic;
		  m->nrl=NL;
		  m->nrh=nrh;
		  m->ncl=NL;
		  m->nch=nch;
		  m->ndl=NL;
		  m->ndh=ndh;


		  m->co=d3tensor(1,ndh,1,nrh,1,nch);


		  return m;


}


/*-----------------------------------------------------------------------*/


void free_d3tensor(double ***t,long ndl, long nrl,  long ncl)

{
ncl=0;
free((FREE_ARG) (t[ndl][nrl]+ndl-NR_END));

free((FREE_ARG) (t[ndl]+nrl-NR_END));
free((FREE_ARG) (t+ndl-NR_END));

}

/*-----------------------------------------------------------------------*/


DOUBLETENSOR* free_doubletensor( DOUBLETENSOR *m)

{


		  if(m==NULL || m->co==NULL){
			  	t_error("This matrix was never allocated");
		}else if(m->isdynamic==1){

			free_d3tensor(m->co,NL,NL,NL);
			m->isdynamic=m->nrl=m->ncl=m->nrh=m->nch=m->ndl=m->ndh=-1;
			free(m);

			return NULL;

		  }else{
			printf("\nWarning::An attemp was made to free a non dynamic tensor\n");

		  }
		  return m;

}



/*---------------------------------------------------------------------------*/
void initialize_doubletensor(DOUBLETENSOR *L,double sign)

{

long i,j,k;

if(L!=NULL){
	if(L->isdynamic==1){
	for(k=1;k<=L->ndh;k++){
		for(i=1;i<=L->nrh;i++){
			for(j=1;j<=L->nch;j++){			
				L->co[k][i][j]=sign;
				
				}
			}
		}
	}else{
		t_error("This tensor was no properly allocated");
	}
}else {
	t_error("A null tensor was addressed");
}
}








