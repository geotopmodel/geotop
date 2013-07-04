#include "turtle.h"
#include "lu.h"
#include "t_datamanipulation.h"
long ind_matrix(long r, long c, LONGVECTOR *Ax, LONGVECTOR *Ai){

	long indbeg, indend, ind;

	if(c==1){
		indbeg=0;
		indend=Ax->co[c];
	}else{
		indbeg=Ax->co[c-1];
		indend=Ax->co[c];
	}

	ind=indbeg;
	do{
		ind++;
	}while(Ai->co[ind]!=r && ind<indend);

	if(Ai->co[ind]==r){
		//printf("%ld %ld ind:%ld\n",r,c,ind);
		return(ind);
	}else{
		//printf("%ld %ld ind:%ld\n",r,c,0);
		return(0);
	}

}




int mat_lu( DOUBLEVECTOR *Ap, LONGVECTOR *Ax, LONGVECTOR *Ai, LONGVECTOR *P ){

	long i, j, k, n=Ax->nh;
	long maxi, tmp;
	double c, c1;
	long p;
	long ind,ind2,ind3;

	n=P->nh;

	for (p=1,i=1; i<=n; i++){
		P->co[i] = i;
	}

	for (k=1; k<=n; k++){
	/*
	* --- partial pivoting ---
	*/
		printf("K:%ld/%ld\n",k,n);

		for (i=k, maxi=k, c=0.0; i<=n; i++){

			ind=ind_matrix(P->co[i],k,Ax,Ai);
			if(ind>0){
				c1=fabs(Ap->co[ind]);
			}else{
				c1=0.0;
			}
			if (c1 > c){
				c = c1;
				maxi = i;
			}
		}

	/*
	*	row exchange, update permutation vector
	*/
		if (k != maxi){
			p++;
			tmp = P->co[k];
			P->co[k] = P->co[maxi];
			P->co[maxi] = tmp;
		}

	/*
	*	suspected singular matrix
	*/

	 	ind=ind_matrix(P->co[k],k,Ax,Ai);
	 	if(ind>0){
	 		if(Ap->co[ind]==0.0){
				t_error("ERR");
				return(-1);
			}
	 	}else{
			t_error("ERR");
	 		return(-1);
	 	}

		for (i=k+1; i<=n; i++){
		/*
		* --- calculate m(i,j) ---
		*/
			ind2=ind_matrix(P->co[i],k,Ax,Ai);
			ind=ind_matrix(P->co[k],k,Ax,Ai);
			if(ind2>0){
				//printf("--> ind:%ld ind2:%ld %f %f\n",ind,ind2,Ap->co[ind],Ap->co[ind2]);
				Ap->co[ind2]=Ap->co[ind2]/Ap->co[ind];
			}

		/*
		* --- elimination ---
		*/
			for (j=k+1; j<=n; j++){
				ind=ind_matrix(P->co[i],j,Ax,Ai);
				ind3=ind_matrix(P->co[k],j,Ax,Ai);
				if(ind>0 && ind2>0 && ind3>0){
					//printf("-> %ld %ld %ld %f %f %f\n",ind,ind2,ind3,Ap->co[ind],Ap->co[ind2],Ap->co[ind3]);
					Ap->co[ind] -= Ap->co[ind2]*Ap->co[ind3];
				}
			}
		}
	}

	return (p);

}




void mat_backsubs1(DOUBLEVECTOR *Ap, LONGVECTOR *Ax, LONGVECTOR *Ai, DOUBLEVECTOR *B, DOUBLEVECTOR *X, LONGVECTOR *P){

	long	ind, i, j, k, n=Ax->nh;
	double	sum;

	for (k=1; k<=n; k++){
		for (i=k+1; i<=n; i++){
			ind=ind_matrix(P->co[i],k,Ax,Ai);
			if(ind>0) B->co[P->co[i]] -= Ap->co[ind] * B->co[P->co[k]];
		}
	}

	ind=ind_matrix(P->co[n],n,Ax,Ai);
	//printf("....n:%ld ind:%ld %f\n",n,ind,Ap->co[ind]);
    X->co[n] = B->co[P->co[n]]/Ap->co[ind];

	for (k=n-1; k>=1; k--){
		sum = 0.0;
		for (j=k+1; j<=n; j++){
			ind=ind_matrix(P->co[k],j,Ax,Ai);
			if(ind>0) sum += Ap->co[ind]*X->co[j];
		}
		ind=ind_matrix(P->co[k],k,Ax,Ai);
		//printf("....k:%ld ind:%ld %f\n",k,ind,Ap->co[ind]);
		X->co[k] = (B->co[P->co[k]]-sum)/Ap->co[ind];
	}

}




void mat_lsolve(DOUBLEVECTOR *Ap, LONGVECTOR *Ax, LONGVECTOR *Ai, DOUBLEVECTOR *B, DOUBLEVECTOR *X)
{
	LONGVECTOR *P;
	long i, n=Ax->nh;

	initialize_doublevector(X, 0.0);

	P=new_longvector(n);
	printf("oki\n");
	i=mat_lu(Ap, Ax, Ai, P );
	printf("oki2\n");
	mat_backsubs1(Ap, Ax, Ai, B, X, P);
	printf("oki3\n");

	free_longvector(P);
}

