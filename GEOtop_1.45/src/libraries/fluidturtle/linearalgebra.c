//#include "turtle.h"
#include "linearalgebra.h"
static float sqrarg;

/*-----------------------------------------------------*/

void convlv(FLOATVECTOR *data,unsigned long n,FLOATVECTOR *respns,unsigned long m,int isign,
             double delta_time,FLOATVECTOR *ans)

{
	unsigned long i,no2;
	float dum,mag2;
	FLOATVECTOR *fft1;
	FLOATVECTOR *fft;

	fft1=new_floatvector(n<<1);
	
	fft=fft1;
	for (i=1;i<=(m-1)/2;i++){
		respns->co[n+1-i]=respns->co[m+1-i];
		}
	for (i=(m+3)/2;i<=n-(m-1)/2;i++){
		respns->co[i]=0.0;
		}
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans->co[i-1]=(fft->co[i-1]*(dum=ans->co[i-1])-fft->co[i]*ans->co[i])/no2*delta_time;
			ans->co[i]=(fft->co[i]*dum+fft->co[i-1]*ans->co[i])/no2*delta_time;
		} else if (isign == -1) {
			if ((mag2=SQR(ans->co[i-1])+SQR(ans->co[i])) == 0.0)
				t_error("Deconvolving at response zero in convlv");
			ans->co[i-1]=(fft->co[i-1]*(dum=ans->co[i-1])+fft->co[i]*ans->co[i])/mag2/no2;
			ans->co[i]=(fft->co[i]*dum-fft->co[i-1]*ans->co[i])/mag2/no2;
		} else t_error("No meaning for isign in convlv");
	}
	ans->co[2]=ans->co[n+1];
	realft(ans,n,-1);

	for(i=1;i<=n;i++){
	   ans->co[i]=ans->co[i];
	}   
	free_floatvector(fft1);
}

/*-----------------------------------------------------*/

void four1(FLOATVECTOR *data,unsigned long nn,int isign)

{
	unsigned long n,mmax,m,j,istep,i;
	float wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data->co[j],data->co[i]);
			SWAP(data->co[j+1],data->co[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data->co[j]-wi*data->co[j+1];
				tempi=wr*data->co[j+1]+wi*data->co[j];
				data->co[j]=data->co[i]-tempr;
				data->co[j+1]=data->co[i+1]-tempi;
				data->co[i] += tempr;
				data->co[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

/*-----------------------------------------------------*/

void ludcmp(SHORTVECTOR *indx, DOUBLEMATRIX *var)

{
short mode;
long i,j,imax,k,n;
double aamax,dum,sum;
DOUBLEVECTOR *vv;
mode=1;
n=var->nrh;
vv=new_doublevector(n);
/*calcolo l'elemento di modulo massimo nella riga i e divido tutti gli elementi 
  per quel numero*/
for(i=1;i<=n;i++){
   aamax=0;
   for(j=1;j<=n;j++){
     if(fabs(var->co[i][j])>aamax)aamax=fabs(var->co[i][j]);
   }
   if(aamax==0){
     printf("aamax e' nullo nella riga (%ld)\n",i);
     aamax=1.0; 
   }
   vv->co[i]=1.0/aamax;
}
/*Questo ciclo va eseguito per tutte le colonne della matrice*/
for(j=1;j<=n;j++){
/*Questo ciclo cambia il valore della matrice modificando secondo*/
   for(i=1;i<=j-1;i++){
      sum=var->co[i][j];
      for(k=1;k<=i-1;k++){
         sum-=var->co[i][k]*var->co[k][j];
      }
      var->co[i][j]=sum;
   }
   aamax=0.0;
   for(i=j;i<=n;i++){
      sum=var->co[i][j];
      for(k=1;k<=j-1;k++){
         sum-=var->co[i][k]*var->co[k][j];
      }
      var->co[i][j]=sum;
      dum=vv->co[i]*fabs(sum);
      if(dum>aamax){
        imax=i;
        aamax=dum;
      }
   }
   if(j!=imax){
     for(i=1;i<=n;i++){
        dum=var->co[imax][i];
        var->co[imax][i]=var->co[j][i];
        var->co[j][i]=dum;
     }
     vv->co[imax]=vv->co[j];
   }
   indx->co[j]=imax;
   if(var->co[j][j]==0){
    var->co[j][j]=TINY;
    //printf("var e' nullo\n");
   }  

   if(j!=n){
     dum=1.0/var->co[j][j];
     for(i=j+1;i<=n;i++){
        var->co[i][j]*=dum;
     }
   }
}
free_doublevector(vv);
}

/*-----------------------------------------------------*/

void lubksb(DOUBLEMATRIX *var, SHORTVECTOR *indx,DOUBLEVECTOR *gam)

{
long i,j,k,ii,n;
double sum;
ii=0;
n=var->nrh;
for(i=1;i<=n;i++){
   k=indx->co[i];
   sum=gam->co[k];
   gam->co[k]=gam->co[i];
   if(ii!=0){
/*calcolo il valore del termine modificato*/
     for(j=ii;j<=i-1;j++){
        sum-=var->co[i][j]*gam->co[j];
     }
   }else if(sum!=0){
     ii=i;
   }
   gam->co[i]=sum;
}
/*Valuto il vettore incognito eseguendo una back substitution*/
for(i=n;i>=1;i--){
   sum=gam->co[i];
   for(j=i+1;j<=n;j++){
      sum-=var->co[i][j]*gam->co[j];
   }
   gam->co[i]=sum/var->co[i][i];
}
}

/*------------------------------------------*/

void realft(FLOATVECTOR *data,unsigned long n,int isign)
{
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	float wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(float) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data->co[i1]+data->co[i3]);
		h1i=c1*(data->co[i2]-data->co[i4]);
		h2r = -c2*(data->co[i2]+data->co[i4]);
		h2i=c2*(data->co[i1]-data->co[i3]);
		data->co[i1]=h1r+wr*h2r-wi*h2i;
		data->co[i2]=h1i+wr*h2i+wi*h2r;
		data->co[i3]=h1r-wr*h2r+wi*h2i;
		data->co[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data->co[1] = (h1r=data->co[1])+data->co[2];
		data->co[2] = h1r-data->co[2];
	} else {
		data->co[1]=c1*((h1r=data->co[1])+data->co[2]);
		data->co[2]=c1*(h1r-data->co[2]);
		four1(data,n>>1,-1);
	}
}

/*-----------------------------------------------------*/

void twofft(FLOATVECTOR *data1,FLOATVECTOR *data2,FLOATVECTOR *fft1,FLOATVECTOR *fft2,unsigned long n)

{

	unsigned long nn3,nn2,jj,j;

	float rep,rem,aip,aim;



	nn3=1+(nn2=2+n+n);

	for (j=1,jj=2;j<=n;j++,jj+=2) {

		fft1->co[jj-1]=data1->co[j];

		fft1->co[jj]=data2->co[j];

	}

	four1(fft1,n,1);

	fft2->co[1]=fft1->co[2];

	fft1->co[2]=fft2->co[2]=0.0;

	for (j=3;j<=n+1;j+=2) {

		rep=0.5*(fft1->co[j]+fft1->co[nn2-j]);

		rem=0.5*(fft1->co[j]-fft1->co[nn2-j]);

		aip=0.5*(fft1->co[j+1]+fft1->co[nn3-j]);

		aim=0.5*(fft1->co[j+1]-fft1->co[nn3-j]);

		fft1->co[j]=rep;

		fft1->co[j+1]=aim;

		fft1->co[nn2-j]=rep;

		fft1->co[nn3-j] = -aim;

		fft2->co[j]=aip;

		fft2->co[j+1] = -rem;

		fft2->co[nn2-j]=aip;

		fft2->co[nn3-j]=rem;

	}

}



/*
___________________________________________________________________________
===========================================================================
                        FUNZIONE       vett_mat
___________________________________________________________________________
*/
/* Questa funzione converte tre vettori in una matrice quadrata tridiagonale.
   Bisogna passare alla funzione i puntatori agli elementi dei vettori
   - d elementi delle diagonale principale;
   - ds elementi della diagonale superiore;
   - di elementi della diagonale inferiore;
   inoltre bisogna passare n, che e' la dimensione della matrice che
   si vuole generare.
   
   La funzione restituisce il puntatore alla matrice.
 */
   
DOUBLEMATRIX *vett_mat (double *d,double *ds,double *di,int n)
{
 int i,j;
 
 DOUBLEMATRIX *A;
 
 A=new_doublematrix(n,n);
 
 for (i=1; i<=n; i++)
 	{ for (j=1; j<=n; j++) A->co[i][j]=0; }
 	
 A->co[1][2]=ds[1];
 A->co[n][n-1]=di[n-1];	
 
 for (i=1; i<=n; i++) A->co[i][i]=d[i];
 for (i=2; i<=n-1; i++) 
     {
      A->co[i][i+1]=ds[i];
      A->co[i][i-1]=di[i-1];	 
 	}	
 	
 return A;
}

/*DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 79
______________________________________________________________________________
==============================================================================
                      FUNZIONE         sprsin
______________________________________________________________________________
*/
/*
  Questa funzione converte una matrice memorizzata nel modo convenzionale
  in un vettore sa[] che contiene solo i valori non nulli della matrice
  e in un vettore ija[] che permette di individuare la posizione originale
  degli elementi di sa[].
  
  La funzione richiede:
  - **a, un puntatore agli elementi della matrice originale;
  - n, dimensione della matrice;
  - thresh, gli elementi della matrice minori di thresh non vengono
            letti;
  - nmax, la lunghezza dei vettori sa[] e ija[].  
*/  
  

void sprsin(double **a, int n, float thresh, long nmax, double sa[],
	long ija[])
{
	int i,j;
	long k;

	for (j=1;j<=n;j++) sa[j]=a[j][j];
	ija[1]=n+2;
	k=n+1;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) {
			if (fabs(a[i][j]) > thresh && i != j) {
				if (++k > nmax) t_error("sprsin: nmax too small");
				sa[k]=a[i][j];
				ija[k]=j;
			}
		}
		ija[i+1]=k+1;
	}
}

/* DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 86-88
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      linbcg                    
_________________________________________________________________________________
*/
/* Questa funzione consente di risolvere di risolvere un sistema lineare del tipo
   A x = b con il metodo iterativo del gradiente coniugato.
   
   Alla funzione devono essere passati i seguenti argomenti:
   - n: dimensione del sistema;
   - sa[] e ija[]: vettori generati dalla funzione sprssin() che memorizzano la 
     matrice;
   - b[]: elementi del vettore dei termini noti;
   - x[]: elementi del vettore soluzione (in ingresso questo vettore deve contenere
     una soluzione di primo tentativo);
   - itol, tol, itmax: parametri gli definiti sopra.
   
   Oltre alla soluzione la funzione calcola anche il numero di iterazione effetuate
   ( iter ) e l'errore commesso ( err ).
*/           


/* commentata per l'errore (dovuto ad una errata definizione di asolve):
error: conflicting types for `asolve'

void linbcg(long n, double b[], double x[], int itol, double tol,
	int itmax, int *iter, double *err, double sa[], long ija[])
{
	void asolve(long n, double b[], double x[], int itrnsp,double sa[], long ija[]);
	void atimes(long n, double x[], double r[], int itrnsp,double sa[], long ija[]);
	double snrm(long n, double sx[], int itol);
	long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	
	DOUBLEVECTOR *P,*PP,*R,*RR,*Z,*ZZ;
	double *p,*pp,*r,*rr,*z,*zz;

	
	P=new_doublevector(n);
	PP=new_doublevector(n);
	R=new_doublevector(n);
	RR=new_doublevector(n);
     Z=new_doublevector(n);
     ZZ=new_doublevector(n);	
	
	p=P->co;
	pp=PP->co;
	r=R->co;
	rr=RR->co;
	z=Z->co;
	zz=ZZ->co;

	*iter=0;
	atimes(n,x,r,0,sa,ija);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	if (itol == 1) {
		bnrm=snrm(n,b,itol);
		asolve(n,r,z,0,sa,ija);
	}
	else if (itol == 2) {
		asolve(n,b,z,0,sa,ija);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0,sa,ija);
	}
	else if (itol == 3 || itol == 4) {
		asolve(n,b,z,0,sa,ija);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0,sa,ija);
		znrm=snrm(n,z,itol);
	} else t_error("illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		asolve(n,rr,zz,1,sa,ija);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(n,p,z,0,sa,ija);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(n,pp,zz,1,sa,ija);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(n,r,z,0,sa,ija);
		if (itol == 1)
			*err=snrm(n,r,itol)/bnrm;
 		else if (itol == 2)
			*err=snrm(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS1*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
	if (*err <= tol) break;
	}

	free_doublevector(P);
	free_doublevector(PP);
	free_doublevector(R);
	free_doublevector(RR);
	free_doublevector(Z);
	free_doublevector(ZZ);
}

//DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 88
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      snrm                    
_________________________________________________________________________________
*/
/* Questa funzione calcola la norma di un vettore con la modalita' specificata 
   dal parametro itol */


double snrm(long n, double sx[], int itol)
{
	long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}










/* DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 88
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      atimes                    
_________________________________________________________________________________
*/
                          


void atimes(long n, double x[], double r[], int itrnsp,double sa[], long ija[])
{
	void dsprsax(double sa[], long ija[], double x[], double b[],
		long n);
	void dsprstx(double sa[], long ija[], double x[], double b[],
		long n);

	if (itrnsp) dsprstx(sa,ija,x,r,n);
	else dsprsax(sa,ija,x,r,n);
}


/*DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 89
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      asolve                   
_________________________________________________________________________________
*/


void asolve(long n, double b[], double x[],double sa[])
{
	long i;

	for(i=1;i<=n;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
}


/* DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 79
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      dsprsax                   
_________________________________________________________________________________
*/
/* Questa funzione moltiplica una matrice, memoririzzata alla maniera di N.R.,
   per un vettore x[]. Il risultato e' un vettore b[].
*/    


void dsprsax(double sa[], long ija[], double x[], double b[], long n)
{
	long i,k;

	if (ija[1] != n+2) t_error("dsprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}


/* DA NUMERICAL RECEPIS IN C. (Second Edition - Cambridge Univ. Press). pag 80
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      dsprstx                  
_________________________________________________________________________________
*/
/* Questa funzione moltiplica la trasposta di una matrice, memoririzzata alla
   maniera di N.R., per un vettore x[]. Il risultato e' un vettore b[].
*/    


void dsprstx(double sa[], long ija[], double x[], double b[], long n)
{
	long i,j,k;
	if (ija[1] != n+2) t_error("mismatched vector and matrix in dsprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}
/*
_________________________________________________________________________________
=================================================================================
                          FUNZIONE      integral                
_________________________________________________________________________________
 Questa funzione calcola l'integrale come sommatoria della funz func 

 */
 float integration(float (*func)(float), float a, float b, int n)
{
     float x=0.0,tnm=0.0,sum=0.0,del=0.0;
     static float s=0.0;
     int it=0,j=0;
     if (n == 1) {
         return (s=0.5*(b-a)*(func(a)+func(b)));
     } else {
         for (it=1,j=1;j<n-1;j++) {
	 it <<= 1;
	 
	 }
	
         tnm=it;
	  del=(b-a)/tnm;                  
         x=a+0.5*del;
         for (sum=0.0,s=0.0,j=1;j<=it;j++,x+=del) {
	 sum += func(x);
	 
	 }
         s=(s+(b-a)*sum/tnm);               
         return s;
     }
}
