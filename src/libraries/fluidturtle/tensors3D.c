#include  "turtle.h"
#include  "tensor3D.h"
/* Note that depth is the first indices and that the indices were pernutated
with respect to NR */

/*-----------------------------------------------------------------------*/


double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	
	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) t_error("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;
	
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) t_error("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) t_error("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/* return pointer to array of pointers to rows */
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


		  m->co=d3tensor(m->ndl,m->ndh,m->nrl,m->nrh,m->ncl,m->nch);


		  return m;


}

DOUBLETENSOR *new_doubletensor0(long ndh,long nrh,long nch)


{
	
	DOUBLETENSOR *m;
	
	m=(DOUBLETENSOR *)malloc(sizeof(DOUBLETENSOR));
	if (!m) t_error("allocation failure in new_doubletensor()");
	m->isdynamic=isDynamic;
	m->nrl=NL;
	m->nrh=nrh;
	m->ncl=NL;
	m->nch=nch;
	m->ndl=0;
	m->ndh=ndh;
	
	
	m->co=d3tensor(m->ndl,m->ndh,m->nrl,m->nrh,m->ncl,m->nch);
	
	
	return m;
	
	
}


/*-----------------------------------------------------------------------*/


void free_d3tensor(double ***t, long nrl, long ncl, long ndl)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}


/*-----------------------------------------------------------------------*/


void free_doubletensor( DOUBLETENSOR *m)

{


		  if(m==NULL || m->co==NULL){
			  	t_error("This matrix was never allocated");
		}else if(m->isdynamic==1){

			free_d3tensor(m->co,m->ndl,m->nrl,m->ncl);
			m->isdynamic=m->nrl=m->ncl=m->nrh=m->nch=m->ndl=m->ndh=-1;
			free(m);

			return;

		  }else{
			printf("\nWarning::An attemp was made to free a non dynamic tensor\n");

		  }


}



/*---------------------------------------------------------------------------*/
void initialize_doubletensor(DOUBLETENSOR *L, double sign)

{

long i,j,k;

if(L!=NULL){
	if(L->isdynamic==1){
		for(k=L->ndl;k<=L->ndh;k++){
			for(i=L->nrl;i<=L->nrh;i++){
				for(j=L->ncl;j<=L->nch;j++){			
					L->co[k][i][j]=sign;
				}
			}
		}
	}else{
		t_error("This tensor was no properly allocated");
	}
}else{
	t_error("A null tensor was addressed");
}
}








