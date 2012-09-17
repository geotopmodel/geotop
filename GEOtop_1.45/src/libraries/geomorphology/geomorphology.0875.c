
#include "geomorphology.0875.h"

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
void sky_view_factor(DOUBLEMATRIX *sky, long N, T_INIT *UV, DOUBLEMATRIX *input, SHORTMATRIX *convess, long novalue)
{
 long i,j,t,m,n,p,q,h,k,r,s; //counter
 double deltateta; //amplitude of the angles in which the horizon is divided
 DOUBLEMATRIX *alfa; //matrices with the angles of the direction
 DOUBLEVECTOR *vv; //vector with the view factor of the current pixel for one of the N parts
 DOUBLEVECTOR *v; //vector with the minimum view factor of the current pixel for one of the N parts
 double vvv; //mean of the sky view for a pixel of the N parts

 if(sky->nrh!=input->nrh) t_error("Sky view factor fatal error, number of rows not consistent");
 if(sky->nch!=input->nch) t_error("Sky view factor fatal error, number of cols not consistent");

 // Computation of the matrix with the angles of the direction
 alfa=new_doublematrix(2*input->nrh-1,2*input->nch-1);
 initialize_doublematrix(alfa,(double)novalue); //initialisation with novalue
 for(i=1;i<=2*input->nrh-1;i++){
      for(j=1;j<=2*input->nch-1;j++){
            if(i<=input->nrh && j<input->nch){
                alfa->co[i][j]=3.0/2.0*Pi+atan(((input->nrh-i)*UV->U->co[1])/
                                    ((input->nch-j)*UV->U->co[1]));
            }
            if(i>input->nrh && j<=input->nch){
                alfa->co[i][j]=Pi+atan(((input->nch-j)*UV->U->co[1])/
                                    ((i-input->nrh)*UV->U->co[1]));
            }
            if(i>=input->nrh && j>input->nch){
                alfa->co[i][j]=Pi/2.0+atan(((i-input->nrh)*UV->U->co[1])/
                                    ((j-input->nch)*UV->U->co[1]));
            }
            if(i<input->nrh && j>=input->nch){
                alfa->co[i][j]=atan(((j-input->nch)*UV->U->co[1])/
                                    ((input->nrh-i)*UV->U->co[1]));
            }
      }
 }

 // Computation of matrix with sky view factor: 
 for(i=1;i<=sky->nrh;i++){
	for(j=1;j<=sky->nch;j++){
		sky->co[i][j]=(double)novalue;
	}
 }

 v=new_doublevector(N);
 vv=new_doublevector(N);
 deltateta=2.0*Pi/N;

 for(i=1;i<=input->nrh;i++){
	for(j=1;j<=input->nch;j++){
		if ((long)input->co[i][j]!=novalue){ //computation only of novalue pixels
		   	for(t=1;t<=N;t++){
		       	v->co[t]=1.0;
		   	}
		   	m=input->nrh-i+1;
		   	n=input->nch-j+1;
		   	p=m+input->nrh-1;
		   	q=n+input->nch-1;
           	for(h=m;h<=p;h++){
		      	for(k=n;k<=q;k++){
		         	for(t=1;t<=N;t++){
		               	if (alfa->co[h][k]>=(t-1)*deltateta && alfa->co[h][k]<t*deltateta){
		                  	r=h-m+1;
		                  	s=k-n+1;
		                  	if (convess->co[r][s]==1 && sqrt(pow((r-i),2)+pow((s-j),2))!=0){
		                       vv->co[t]=1-sin(atan((input->co[r][s]-input->co[i][j])
		                         	          /(sqrt(pow((r-i),2)+pow((s-j),2))*UV->U->co[1])));
                             	if (vv->co[t]<v->co[t]){
                                   v->co[t]=vv->co[t];
                             	}
                          	}
                          	break;
                       	}
                 	}
              	}
           	}
           	vvv=0.0;
           	for(t=1;t<=N;t++){
               	vvv=vvv+v->co[t];
           	}
           	sky->co[i][j]=(1.0/N*vvv);
         }
   	}
   	printf("Percentage of the calculation of the sky view factor matrix: %5.2f%%\n",
   	       100.0*(double)i/(double)sky->nrh);
}

free_doublematrix(alfa);
free_doublevector(v);
free_doublevector(vv);
	
}

//***************************************************************************


void pits_filler_0875(DOUBLEMATRIX *Z0,SHORTMATRIX *land_use)
//! Subroutine to fill pits in a DEM
{
long i,j,k,ii,jj,nx,ny,N1=0,N2=0,mod=0,ISPIT,m,inb,jnb,clustdim,dimlistaor,pp,ne;
long imin,jmin,Nsfioro,iP,jP,Ncand,indP,indQ,iQ,jQ,contr;
long dim=100000; /*dimensione consentita dei vettori contenti pits*/
double zclust,zminInt,zminInt_1,rndNum,incr;
long data [9][2]={ {0,0},
                   {0,1},
                   {-1,1},
                   {-1,0},
                   {-1,-1},
                   {0,-1},
                   {1,-1},
                   {1,0},
                   {1,1}, };
LONGVECTOR *mod_i,*mod_j,*is,*js,*iclust,*jclust,*ior,*jor,*indcand;
DOUBLEVECTOR *nearel;

/*CALCOLO ORA QUANTI PITS INIZIALI SONO PRESENTI NELLA MIA MATRICE*/
nx=Z0->nrh;
ny=Z0->nch;
ne=nx*ny;
mod_i=new_longvector(ne);
mod_j=new_longvector(ne);
initialize_longvector(mod_i,0);
initialize_longvector(mod_j,0);
N2=0;
do{
N1=0;
is=new_longvector(dim);
js=new_longvector(dim);
ISPIT=0;
for(i=1;i<=nx;i++){
   for(j=1;j<=ny;j++){
     if(Z0->co[i][j]>=10000){
       Z0->co[i][j]=Z0->co[1][1];/*Si suppone che il primo pixel sia un NOVALUE*/
     }
   }
}
for(i=2;i<=nx-1;i++){
for(j=2;j<=ny-1;j++){
   /*Escludo dal controllo i punti del contorno (con Z0=NOVALUE o =10000) e le uscite
     (con elevazione posta a 10^-3)*/
   if (Z0->co[i][j]>0.0011 && Z0->co[i][j]<10000){
      nearel=new_doublevector(9);
      for(k=1;k<=8;k++){
         ii=i+data[k][0];
         jj=j+data[k][1];
         nearel->co[k]=Z0->co[ii][jj];
      }
      nearel->co[9]=Z0->co[i][j];
      if(nearel->co[1]>=nearel->co[9] &&
      nearel->co[2]>=nearel->co[9] &&
      nearel->co[3]>=nearel->co[9] &&
      nearel->co[4]>=nearel->co[9] &&
      nearel->co[5]>=nearel->co[9] &&
      nearel->co[6]>=nearel->co[9] &&
      nearel->co[7]>=nearel->co[9] &&
      nearel->co[8]>=nearel->co[9] &&
      land_use->co[i][j]!=11){
         ISPIT=1;
      }else{
          ISPIT=0;
      }
      free_doublevector(nearel);
      if (ISPIT==1){
         N1+=1;
         if (N1>dim){
            printf("\nAumentare il valore della variabile interna dim nella ");
            printf("\nsubroutine pits_filler_0875.\n");
         }
         is->co[N1]=i;
         js->co[N1]=j;
      }
   }
}
}
N2+=1;
printf("Il numero di N (%ld)(%ld)\n",N1,N2);
if(N1!=0){
for(m=1;m<=N1;m++){
   iclust=new_longvector(dim);
   jclust=new_longvector(dim);
   ior=new_longvector(dim);
   jor=new_longvector(dim);
   indcand=new_longvector(10);
   /*Definisco il cluster iniziale comprendente solo il pit*/
   zclust=Z0->co[is->co[m]][js->co[m]];
   iclust->co[1]=is->co[m];
   jclust->co[1]=js->co[m];
   /*Valuto se il pit e' effettivamente ancora un pit*/
   nearel=new_doublevector(9);
   for(k=1;k<=8;k++){
      i=is->co[m]+data[k][0];
      j=js->co[m]+data[k][1];
      nearel->co[k]=Z0->co[i][j];
   }
   nearel->co[9]=Z0->co[is->co[m]][js->co[m]];
   if (nearel->co[1]+0.00001>=nearel->co[9] &&
   nearel->co[2]+0.00001>=nearel->co[9] &&
   nearel->co[3]+0.00001>=nearel->co[9] &&
   nearel->co[4]+0.00001>=nearel->co[9] &&
   nearel->co[5]+0.00001>=nearel->co[9] &&
   nearel->co[6]+0.00001>=nearel->co[9] &&
   nearel->co[7]+0.00001>=nearel->co[9] &&
   nearel->co[8]+0.00001>=nearel->co[9] &&
   land_use->co[is->co[m]][js->co[m]]!=11){
      ISPIT=1;
   }else{
      ISPIT=0;
   }
   free_doublevector(nearel);
   if (ISPIT==1){
      clustdim=1;
      /*Valuto il punto più basso dell'intorno del cluster.*/
       zminInt=10000;
       for(k=0;k<=7;k++){
           contr=0;
           inb=iclust->co[1]+data[k][0];
           if(inb==0 || inb==(nx+1)) contr=1;
           jnb=jclust->co[1]+data[k][1];
           if(jnb==0 || jnb==(ny+1)) contr=1;
           /*Controllo se il vicino del punto del cluster e' esso stesso un punto del cluster
             (in tal caso passo all'altro vicino)*/
           for(j=1;j<=clustdim;j++){
              if (inb==iclust->co[j] && jnb==jclust->co[j]){
                 contr=1;
              }
           }
           if (Z0->co[inb][jnb]<0.0009) contr=1;
           if (contr==0){
              if (Z0->co[inb][jnb] < zminInt){
                 zminInt=Z0->co[inb][jnb];
                 imin=inb;
                 jmin=jnb;
              }
           }
        }
        zminInt_1=zminInt;
        clustdim+=1;
        iclust->co[clustdim]=imin;
        jclust->co[clustdim]=jmin;
        Nsfioro=clustdim;
        /*Procedo all'aggiustamento delle direzioni di drenaggio:
        1. Prendo il primo punto Nsfioro.*/
        dimlistaor=1;
        ior->co[1]=iclust->co[Nsfioro];
        jor->co[1]=jclust->co[Nsfioro];
        clustdim-=1;
        /* Elimino Nsfioro dalla lista dei pixel appartenenti al cluster
           (se non si dovesse trovare nell'ultima posizione vanno spostati tutti
           gli elementi a destra di Nsfioro di un posto a sinistra)*/
        if(Nsfioro < (clustdim+1)){
           for(i=Nsfioro;i<=clustdim;i++){
              j=i+1;
              iclust->co[i]=iclust->co[j];
              jclust->co[i]=jclust->co[j];
           }
        }
        /*2. Scelgo P come punto iniziale il punto 1 della listaor. Il do e' necessario
             in quanto bisogna eliminare i punti origine che non presentano candidati
             Ncand==0(punti intorno al pixel appartenenti al cluster)*/
        do{
           do{
              pp=1313;
              rndNum=ran1(&pp);
              if (rndNum==1.){
                 rndNum=1.0-1.E-3;
              }
              indP=floor(rndNum*dimlistaor)+1;
              iP=ior->co[indP];
              jP=jor->co[indP];
        /*3. Creo la lista dei punti all'intorno di P appartenenti al cluster;
          uno di questi punti (definiti Q) sarà la destinazione della crescita
          di Eden*/
              Ncand=0;
              for(k=1;k<=8;k++){
                 inb=iP+data[k][0];
                 jnb=jP+data[k][1];
                 for(i=1;i<=clustdim;i++){
                    if (inb==iclust->co[i] && jnb==jclust->co[i]){
                       Ncand+=1;
                       indcand->co[Ncand]=i;
                    }
                 }
              }
              /*Se Ncand==0 elimino P dalla lista origine e torno al punto 2*/
              if (Ncand==0){
                 dimlistaor-=1;
                 if (indP < (dimlistaor+1)){
                    for(i=indP;i<=dimlistaor;i++){
                       j=i+1;
                       ior->co[i]=ior->co[j];
                       jor->co[i]=jor->co[j];
                    }
                 }
              }
           }while(Ncand==0);
           /* Prendo il primo punto Q della lista dei candidati ed incremento
              la quota di quel punto */
           pp=1313;
           rndNum=ran1(&pp);
           if (rndNum==1.){
              rndNum=1.0-1.E-3;
           }
           indQ=floor(rndNum*Ncand)+1;
           i=indcand->co[indQ];
           iQ=iclust->co[i];
           jQ=jclust->co[i];
           indQ=floor(rndNum*Ncand)+1;
           i=indcand->co[indQ];
           iQ=iclust->co[i];
           jQ=jclust->co[i];
           incr=0.001;
           if(N2>20)incr+=0.1;
           Z0->co[iQ][jQ]=Z0->co[iP][jP]+incr;
           for(jj=1;jj<=mod;jj++){
              if (iQ==mod_i->co[jj] && jQ==mod_j->co[jj]) contr=1;
           }
           if (contr!=1){
              mod+=1;
              mod_i->co[mod]=iQ;
              mod_j->co[mod]=jQ;
           }
           contr=0;
           /*Aggiungo Q nella listaor*/
           dimlistaor+=1;
           ior->co[dimlistaor]=iQ;
           jor->co[dimlistaor]=jQ;
           /*Elimino Q dalla lista del cluster*/
           clustdim-=1;
           if (indcand->co[indQ]<(clustdim+1)){
              for(i=indcand->co[indQ];i<=clustdim;i++){
                 j=i+1;
                 iclust->co[i]=iclust->co[j];
                 jclust->co[i]=jclust->co[j];
              }
           }
           /*3. Ora se il clustdim!=0( non ho eliminato tutto il cluster) passo a quello sucessivo
                ritorno al punto1 altrimenti passo al cluster sucessivo*/
        }while(clustdim!=0);
     }
     free_longvector(iclust);
     free_longvector(jclust);
     free_longvector(ior);
     free_longvector(jor);
     free_longvector(indcand);
}
}
free_longvector(is);
free_longvector(js);
}while(N1!=0);
free_longvector(mod_i);
free_longvector(mod_j);
}


//***************************************************************************









//Presa da geomorphology099 e modificato der_min
void nablaquadro_mask(DOUBLEMATRIX *Z0,SHORTMATRIX *curv,DOUBLEVECTOR *U,DOUBLEVECTOR *V)

{

short y;
long i,j,h,rows,cols;
double grid[9],z[9],derivate2;
double der_min=0.00001; /*limite per la limite per la planarita'*/

short v[13][2] = {           { 0, 0},
                             { 0, 1},
                             {-1, 1},
                             {-1, 0},
                             {-1,-1},
                             { 0,-1},
                             { 1,-1},
                             { 1, 0},
                             { 1, 1},
                             { 0, 0},
                             { 0, 0},
                             { 0, 0},
                             { 0, 0}             };

grid[0]=0;
grid[1]=grid[5]=U->co[1];
grid[3]=grid[7]=U->co[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

rows=Z0->nrh;
cols=Z0->nch;

for(i=2;i<=rows-1;i++){
   	for(j=2;j<=cols-1;j++){
   		z[0]=Z0->co[i][j];
   		if(z[0]!=V->co[2]){
       		y=1;
        	for(h=1;h<=8;h++){
        		z[h]=Z0->co[i+v[h][0]][j+v[h][1]];
				if(z[h]==V->co[2]){
					y=0;
		        	break;
        		}
       		}
        	if(y==0){
        		curv->co[i][j]=1;
        	}else{
		    	derivate2=0.5*((z[1]+z[5]-2*z[0])/(grid[1]*grid[1])+ (z[3]+z[7]-2*z[0])/(grid[3]*grid[3]));
		    	derivate2=derivate2+0.5*((z[2]+z[4]+z[6]+z[8]-4*z[0])/(grid[6]*grid[6]));

		    	if(fabs(derivate2)<=der_min || derivate2>der_min){  //plane or concave
		    		curv->co[i][j]=0;
		    	}else{
		    		curv->co[i][j]=1;	//convex
		    	}
		 	}
    	}
	}
}
}










//***************************************************************************









/* Calculation of the aspect for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - azimuth   matrix with the aspect
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   program aspetto of Pegoretti; this subroutine is more indipendent and needs
   less input                                                                   */
void aspect0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *azimuth)
{
 long i,j;
 double delta;
 DOUBLEMATRIX *a;/*contiene, per ogni punto del bacino, il coefficiente angolare della retta
                 interpolante il punto a monte, quello a valle e il punto stesso*/
 DOUBLEMATRIX *b;/*contiene, per ogni punto del bacino,il coefficiente angolare della retta
                 interpolante il punto a destra, quello a sinistra e il punto stesso*/
 DOUBLEMATRIX *slope;/*contiente la pendenza media di ogni pixel del bacino*/

 a=new_doublematrix(Z0->nrh,Z0->nch);
 b=new_doublematrix(Z0->nrh,Z0->nch);
 slope=new_doublematrix(Z0->nrh,Z0->nch);
 
 initialize_doublematrix(a,V->co[2]);

 for(i=2;i<=Z0->nrh-1;i++){
      for(j=2;j<=Z0->nch-1;j++){
           if(Z0->co[i][j]!=V->co[2]){
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i+1][j])/(2*U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i][j]-Z0->co[i+1][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=V->co[2];
               }
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j+1])/(2*U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j]-Z0->co[i][j+1])/(U->co[1]));
               }
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j])/(U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=0;
               }
          }
     }
 }

 /*La matrice "slope" contiene la pendenza di ogni pixel:*/
 for(i=1;i<=Z0->nrh;i++){
      for(j=1;j<=Z0->nch;j++){
           if(a->co[i][j]!=V->co[2]){
              	slope->co[i][j]=acos(cos(fabs(a->co[i][j]))*cos(fabs(b->co[i][j])));
           }else{
           		slope->co[i][j]=V->co[2];
           }
      }
 }

 /*La matrice "azimuth" contiene la direzione rispetto al nord di ogni pixel:*/
 for(i=1;i<=Z0->nrh;i++){
  	for(j=1;j<=Z0->nch;j++){
     	if(a->co[i][j]!=V->co[2]){
         	if(a->co[i][j]<0 && b->co[i][j]>0){
           		delta=acos(sin(fabs(a->co[i][j]))*cos(fabs(b->co[i][j]))/
                	(sqrt((double)1-pow(cos(a->co[i][j]),(double)2)*pow(cos(b->co[i][j]),(double)2))));
            	azimuth->co[i][j]=delta;
    		}
 			if(a->co[i][j]>0 && b->co[i][j]>0){
   				delta=acos(sin(fabs(a->co[i][j]))*cos(fabs(b->co[i][j]))/
    				(sqrt((double)1-pow(cos(a->co[i][j]),(double)2)*pow(cos(b->co[i][j]),(double)2))));
          		azimuth->co[i][j]=Pi-delta;
         	}
      		if(a->co[i][j]>0 && b->co[i][j]<0){
          		delta=acos(sin(fabs(a->co[i][j]))*cos(fabs(b->co[i][j]))/
                    (sqrt((double)1-pow(cos(a->co[i][j]),(double)2)*pow(cos(b->co[i][j]),(double)2))));
          		azimuth->co[i][j]=Pi+delta;
        	}
      		if(a->co[i][j]<0 && b->co[i][j]<0){
          		delta=acos(sin(fabs(a->co[i][j]))*cos(fabs(b->co[i][j]))/
     				(sqrt((double)1-pow(cos(a->co[i][j]),(double)2)*pow(cos(b->co[i][j]),(double)2))));
          		azimuth->co[i][j]=2*Pi-delta;
        	}
      		if(a->co[i][j]==0 && b->co[i][j]>0){
          		azimuth->co[i][j]=Pi/2;
        	}
      		if(a->co[i][j]==0 && b->co[i][j]<0){
          		azimuth->co[i][j]=Pi*3./2.;
        	}
      		if(a->co[i][j]>0 && b->co[i][j]==0){
          		azimuth->co[i][j]=Pi;
        	}
      		if(a->co[i][j]<0 && b->co[i][j]==0){
          		azimuth->co[i][j]=0;
        	}
        	if(a->co[i][j]==0 && b->co[i][j]==0){
          		azimuth->co[i][j]=0;
        	}
    	}else{
           	azimuth->co[i][j]=V->co[2];
        }
    }
 }

 free_doublematrix(a);
 free_doublematrix(b);
 free_doublematrix(slope);

}













//***************************************************************************




















/* Calculation of the mean slope for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - slope    matrix with the mean slope [radiants]
   Subroutine created by Davide Tamanini (June 2003)                                */
void slopes0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *slope)
{
 long i,j;
 DOUBLEMATRIX *a;/*contiene, per ogni punto del bacino, il coefficiente angolare della retta
                 interpolante il punto a monte, quello a valle e il punto stesso*/
 DOUBLEMATRIX *b;/*contiene, per ogni punto del bacino,il coefficiente angolare della retta
                 interpolante il punto a destra, quello a sinistra e il punto stesso*/

 a=new_doublematrix(Z0->nrh,Z0->nch);
 b=new_doublematrix(Z0->nrh,Z0->nch);
 
 for(i=2;i<=Z0->nrh-1;i++){
      for(j=2;j<=Z0->nch-1;j++){
           if(Z0->co[i][j]!=V->co[2] ){
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i+1][j])/(2*U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i][j]-Z0->co[i+1][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=0.0;
               }
			   
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j+1])/(2*U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j]-Z0->co[i][j+1])/(U->co[1]));
               }
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j])/(U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=0.0;
               }
          }
     }
 }

 /*La matrice "slope" contiene la pendenza di ogni pixel:*/
 for(i=1;i<=Z0->nrh;i++){
      for(j=1;j<=Z0->nch;j++){
			slope->co[i][j]=acos(cos(fabs(a->co[i][j]))*cos(fabs(b->co[i][j])));
      }
 }

 free_doublematrix(a);
 free_doublematrix(b);

}






//***************************************************************************








//presa uguale da geomorphology099
void nablaquadro(DOUBLEMATRIX *Z0,DOUBLEMATRIX *nabla,DOUBLEVECTOR *U,DOUBLEVECTOR *V)
{
short y,mode;
long i,j,h,n;
double grid[9], z[9],nablaT;
SHORTMATRIX *segn;
short v[9][2] = {
                             { 0, 0},
                             { 0, 1},
                             {-1, 1},
                             {-1, 0},
                             {-1,-1},
                             { 0,-1},
                             { 1,-1},
                             { 1, 0},
                             { 1, 1},
                            };
mode=0;
grid[0]=0;
grid[1]=grid[5]=U->co[1];
grid[3]=grid[7]=U->co[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);
segn=new_shortmatrix(Z0->nrh,Z0->nch);
initialize_shortmatrix(segn,0);
for(i=2;i<=Z0->nrh-1;i++){
       	for(j=2;j<=Z0->nch-1;j++){
           	z[0]=Z0->co[i][j];
            if(z[0]!=V->co[2]) {
	        	y=1;
	        	for(h=1;h<=8;h++){
	        		z[h]=Z0->co[i+v[h][0]][j+v[h][1]];
					if(z[h]==V->co[2]){
				        y=0;
					    segn->co[i][j]=1;
	                	break;
	        	    }
	          	}
	        	if(y==0){
	         		nabla->co[i][j]=V->co[2];
	           	}else{
					nabla->co[i][j]=0.5*((z[1]+z[5]-2*z[0])/(grid[1]*grid[1])+ (z[3]+z[7]-2*z[0])/(grid[3]*grid[3]));
					nabla->co[i][j]+=0.5*((z[2]+z[4]+z[6]+z[8]-4*z[0])/(grid[6]*grid[6]));
               	}
			}else{
				nabla->co[i][j]=V->co[2];
			}
     	}
	}
	for(i=2;i<=Z0->nrh-1;i++){
        for(j=2;j<=Z0->nch-1;j++){
           	if(segn->co[i][j]==1){
            	n=0;
            	nablaT=0;
            	y=0;
	        	for(h=1;h<=8;h++){
        		    z[h]=Z0->co[i+v[h][0]][j+v[h][1]];
        		    y=0;
			        if(z[h]==V->co[2] || nabla->co[i+v[h][0]][j+v[h][1]]==V->co[2]) y=1;
			        if(y==0){
                        n+=1;
                        nablaT+=nabla->co[i+v[h][0]][j+v[h][1]];
                    }
	        	}
				if(n==0)n=1;
            	nabla->co[i][j]=nablaT/(double)n;
           	}
      	}
	}
free_shortmatrix(segn);
}




























//***************************************************************************




/*Computation of the thickness of hydrological active soil (bedrock depth) using
  the method of Dietrich:
   Input:  - soil_parameters  par of soil (in particular that of Dietrich)
           - UV               struct with vector of pixel dimensions and novalue
           - Z0               matrix with elevation (DTM)
   Output: - h                matrix with the depth of bedrock [mm]
   Subroutine created by Davide Tamanini (June 2003) on the basis of the
   subroutine profondita of Pegoretti; this subroutine is more indipendent
   and needs less input; it was modified by Davide Tamanini (April 2004)
   to eliminate the matrix of dranaige direction as input.                      */
void soil_depth(double min_h, double h_crit, double P, double M, double prs, double k , T_INIT *UV, DOUBLEMATRIX *Z0, DOUBLEMATRIX *h)
{
 long i,j,NOVALUES;
 double por,Dcr;
 DOUBLEMATRIX *D_z_sec;/*matrix with the value of second derivative of elevation*/

 por=1-(1/prs);
 Dcr=P*prs/k;
 NOVALUES=(long)UV->V->co[2];

 /*Computation of matrix with the value of second derivative of elevation:*/
 D_z_sec=new_doublematrix(Z0->nrh,Z0->nch);
 nablaquadro(Z0,D_z_sec,UV->U,UV->V);

 /*Computation of bedrock depth:*/
 for(i=2;i<=D_z_sec->nrh-1;i++){
   	for(j=2;j<=D_z_sec->nch-1;j++){
   	   if (D_z_sec->co[i][j]!=UV->V->co[2]){
   		  /* con curv<0 h assegnata con Dietrich */
      	  if (D_z_sec->co[i][j]!=NOVALUES && D_z_sec->co[i][j]<0){
             h->co[i][j]=1000.0/M*(log(Dcr)-log(-(D_z_sec->co[i][j])));
             if(h->co[i][j]<=min_h) h->co[i][j]=min_h;
             if(h->co[i][j]>h_crit) h->co[i][j]=h_crit;
      	  }
      	  if (D_z_sec->co[i][j]!=NOVALUES && D_z_sec->co[i][j]>=0){
        	 h->co[i][j]=h_crit;
      	  }
       }else{
          h->co[i][j]=UV->V->co[2];
       }
   	}
 }

 free_doublematrix(D_z_sec);

}



//***************************************************************************













/* Computation of the Drainage direction following this convention:
   0 for no draining pixel; from 1 to 8 starting from east in counterclockwise toward;
   9 for novalue pixel; 10 for outlet pixel; 11 for lake and sea; 12 pit pixel
   Inputs   - Z0               matrix with elevation (DTM)
            - land_use         matrix with land use
            - UV               struct with vector of pixel dimensions and novalue
   Outputs:	- directions       matrix with Drainage Direction                         */
void DrainageDirections0875(DOUBLEMATRIX *elevations,SHORTMATRIX *land_use,T_INIT *UV,
                            SHORTMATRIX *directions)
{
/* Codici di directions:
0: non drena (dopo il passaggio di "select_hillslopes" sono 0 anche i pixel canale ?)
1: da 1 a 8 drena nei pixels adiacenti (da E antiorario)
2:
3:
4:
5:
6:
7:
8:
9:	novalue
10:	punto di uscita
11: lago (segna=1)
12: punto di pit
*/
 long i,j,k;
 static long dir[11][3] = {	{ 0, 0, 0},
							{ 0, 1, 1},
 							{-1, 1, 2},
 							{-1, 0, 3},
 							{-1,-1, 4},
 							{ 0,-1, 5},
 							{ 1,-1, 6},
 							{ 1, 0, 7},
 							{ 1, 1, 8},
							{ 0, 0, 9},
 							{ 0, 0, 10} };
 double grid[9];
 double df=0.0,steepestdescent=0.0;
 short cont=0;

 grid[0]=1;
 grid[1]=grid[5]=UV->U->co[1];
 grid[3]=grid[7]=UV->U->co[2];
 grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

 initialize_shortmatrix(directions,9);
 for(i=2;i<=elevations->nrh-1;i++){
	for(j=2;j<=elevations->nch-1; j++){
		if(elevations->co[i][j]!=UV->V->co[2] && (land_use->co[i][j]!=11 ||
		                                                    land_use->co[i][j]!=12)){
			steepestdescent=0.0;
			for(k=1;k<=8;k++){
				if(elevations->co[i+dir[k][0]][j+dir[k][1]]!=UV->V->co[2]){
					df=(elevations->co[i][j]-elevations->co[i+dir[k][0]][j+dir[k][1]])/grid[k];
					if(df> steepestdescent){
						steepestdescent=df;
						directions->co[i][j]=k;
					}
				}
			}
			if(directions->co[i][j]==9){
			/*Call of is_ontheborder(elevations,UV->V,i,j):*/
				if (is_ontheborder(elevations,UV->V,i,j)){
					directions->co[i][j]=10;
					printf("Outlet at position %ld %ld\n", i,j);
					cont++;
			    }else{
					printf("1Warning::a pit was found at position %ld %ld\n",i,j);
					directions->co[i][j]=12;
			    }
		   	}
		}
		if(land_use->co[i][j]==11 || land_use->co[i][j]==12) directions->co[i][j]=11;
 	}
}
if(cont>1) printf("Warning::this matrix has multiple outlets \n");
}






//***************************************************************************









//presa uguale da geomorphology099 cambiato il novalue
void gradients(DOUBLEMATRIX *Z0,SHORTMATRIX *directions,DOUBLEMATRIX *dr_gradient,
               T_INIT *UV)

{
 long i,j,rows,cols,k,point[2];
 short     dir[12][3] =      {	{ 0, 0, 0},
								{ 0, 1, 1},
								{-1, 1, 2},
								{-1, 0, 3},
								{-1,-1, 4},
								{ 0,-1, 5},
								{ 1,-1, 6},
								{ 1, 0, 7},
								{ 1, 1, 8},
								{ 0, 0, 9},
								{ 0, 0, 10},
								{ 0, 0, 11}	 };
 double grid[11];

 grid[0]=0;
 grid[1]=grid[5]=UV->U->co[1];
 grid[3]=grid[7]=UV->U->co[2];
 grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);
 grid[10]=UV->U->co[1];

 rows=Z0->nrh;
 cols=Z0->nch;
 for(i=2;i<=rows-1;i++){
	for(j=2;j<=cols-1;j++){
		if  (directions->co[i][j]>=1 && directions->co[i][j]<=8){
			k=directions->co[i][j];
			point[0]=i+dir[k][0];
			point[1]=j+dir[k][1];
			dr_gradient->co[i][j]=(Z0->co[i][j]-Z0->co[point[0]][point[1]])/grid[k];
		}
	}
 }

}



//***************************************************************************







/* Calculation of area considering the slope for each pixels:
   Input:  - Z0        matrix with elevation (DTM)
           - U         vector with pixel dimensions and coordinate of the basin
           - V         vector with the novalues
   Output: - area      matrix with the mean slope
   Subroutine created by Davide Tamanini (June 2003)                                */
void area0875(DOUBLEMATRIX *Z0,DOUBLEVECTOR *U,DOUBLEVECTOR *V,DOUBLEMATRIX *area)
{
 long i,j;
 DOUBLEMATRIX *a;/*contiene, per ogni punto del bacino, il coefficiente angolare della retta
                 interpolante il punto a monte, quello a valle e il punto stesso*/
 DOUBLEMATRIX *b;/*contiene, per ogni punto del bacino,il coefficiente angolare della retta
                 interpolante il punto a destra, quello a sinistra e il punto stesso*/

 a=new_doublematrix(Z0->nrh,Z0->nch);
 b=new_doublematrix(Z0->nrh,Z0->nch);
 
 initialize_doublematrix(a,V->co[2]);

 for(i=2;i<=Z0->nrh-1;i++){
      for(j=2;j<=Z0->nch-1;j++){
           if(Z0->co[i][j]!=V->co[2]){
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i+1][j])/(2*U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]!=V->co[2]){
                  a->co[i][j]=atan((Z0->co[i][j]-Z0->co[i+1][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]!=V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=atan((Z0->co[i-1][j]-Z0->co[i][j])/(U->co[2]));
               }
               if(Z0->co[i-1][j]==V->co[2] && Z0->co[i+1][j]==V->co[2]){
                  a->co[i][j]=V->co[2];
               }
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j+1])/(2*U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]!=V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j]-Z0->co[i][j+1])/(U->co[1]));
               }
               if(Z0->co[i][j-1]!=V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=atan((Z0->co[i][j-1]-Z0->co[i][j])/(U->co[1]));
               }
               if(Z0->co[i][j-1]==V->co[2] && Z0->co[i][j+1]==V->co[2]){
                  b->co[i][j]=0;
               }
          }
     }
 }

 /*La matrice "slope" contiene la pendenza di ogni pixel:*/
 for(i=1;i<=Z0->nrh;i++){
      for(j=1;j<=Z0->nch;j++){
           if (a->co[i][j]!=V->co[2]){
              area->co[i][j]=(U->co[2]/fabs(cos(a->co[i][j])))*(U->co[1]
                                  /fabs(cos(b->co[i][j])));
           }else{
           	  area->co[i][j]=V->co[2];
           }
      }
 }

 free_doublematrix(a);
 free_doublematrix(b);

}


//***************************************************************************




/* Calculation of channel network: this subroutine is a modification of select_hillslopes
   created by Davide Tamanini (August 2003)                                                */
void select_hillslopes_mod(LONGMATRIX *ca,DOUBLEMATRIX *dr_gradient,DOUBLEMATRIX *curv,
                           SHORTMATRIX *m,double threshold,DOUBLEVECTOR *pixelsize)

{
long i,j,flw[2],counter=0,hillslopecounter=0,max;
double laplacian,slope,valore,area;
SHORTMATRIX *flow;
short mode;

max=0;
mode=1;
printf("threshold(%f)\n",threshold);
flow=new_shortmatrix(m->nrh,m->nch);
copy_shortmatrix(m,flow);
	for(i=2;i<=m->nrh-1;i++){
		for(j=2;j<=m->nch-1;j++){
            if (m->co[i][j]<9){
               laplacian=curv->co[i][j];
               slope=dr_gradient->co[i][j];
               area=(double)ca->co[i][j]*pixelsize->co[1]*pixelsize->co[2];
               valore=sqrt((double)area)*(double)slope;
			   if (valore < threshold || laplacian < 0){
			      m->co[i][j]=0;
			   }else{
			      m->co[i][j]=-m->co[i][j];
			   }
			}
		}
	}

	for(i=1;i<=m->nrh;i++){
		for(j=1;j<=m->nch;j++){
			if(m->co[i][j] < 0){
			    counter++;
				flw[0]=i;
				flw[1]=j;
				m->co[i][j]=10;
				go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);

				while( flow->co[flw[0]][flw[1]] < 9 && m->co[flw[0]][flw[1]]==0){
			        counter++;
			   		m->co[flw[0]][flw[1]]=10;
			    	/*printf("%d %d %d\n", m->co[flw[0]][flw[1]],flw[0],flw[1]);*/
					go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);
					/* scanf("%c",&ch); */
				}
			}
		}
	}

printf("Channel pixels are: %ld\n",counter);
	for(i=1;i<=m->nrh;i++){
		for(j=1;j<=m->nch;j++){

			if(m->co[i][j] == 0){
			    hillslopecounter++;
				flw[0]=i;
				flw[1]=j;
				m->co[i][j]=flow->co[i][j];
				go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);

				while(m->co[flw[0]][flw[1]]==0){
			        hillslopecounter++;
			   		m->co[flw[0]][flw[1]]=flow->co[flw[0]][flw[1]];
					go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);
				}
			}
		}
	}

free_shortmatrix(flow);
printf("Hillslope pixels are: %ld\n",hillslopecounter);

}





//***************************************************************************




/*Computation of the trasversal slope of the channel pixel (created by Davide 30-8-2003):*/
void channel_lateral_slope(DOUBLEMATRIX *Z0,SHORTMATRIX *DD,T_INIT *UV,DOUBLEMATRIX *i_ch)
{
 long i,j,rows,cols,k,point[2];
 short     dir[12][3] =      {	{ 0, 0, 0},
								{ 0, 1, 1},
								{-1, 1, 2},
								{-1, 0, 3},
								{-1,-1, 4},
								{ 0,-1, 5},
								{ 1,-1, 6},
								{ 1, 0, 7},
								{ 1, 1, 8},
								{ 0, 0, 9},
								{ 0, 0, 10},
								{ 0, 0, 11}	 };
 double grid[11];
 grid[0]=0;
 grid[1]=grid[5]=UV->U->co[1];
 grid[3]=grid[7]=UV->U->co[2];
 grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);
 grid[10]=UV->U->co[1];
 rows=Z0->nrh;
 cols=Z0->nch;

 for(i=2;i<=rows-1;i++){
	for(j=2;j<=cols-1;j++){
	   i_ch->co[i][j]=0.0001;
		if  (DD->co[i][j]>0 && DD->co[i][j]<9){
			k=DD->co[i][j]+2;
			if (k>8) k-=8;
			point[0]=i+dir[k][0];
			point[1]=j+dir[k][1];
			if (Z0->co[point[0]][point[1]]!=(double)UV->V->co[2]){
			   i_ch->co[i][j]=fabs((Z0->co[i][j]-Z0->co[point[0]][point[1]])
			                       /(grid[k]*2.0));
			}else{
			   i_ch->co[i][j]=-1.0;
			}
			k=DD->co[i][j]-2;
			if (k<1) k+=8;
			point[0]=i+dir[k][0];
			point[1]=j+dir[k][1];
			if (Z0->co[point[0]][point[1]]!=(double)UV->V->co[2]){
			   if (i_ch->co[i][j]<0.0){
			      i_ch->co[i][j]=fabs((Z0->co[i][j]-Z0->co[point[0]][point[1]])
			                          /grid[k]);
			   }else{
			   i_ch->co[i][j]+=fabs((Z0->co[i][j]-Z0->co[point[0]][point[1]])/
			                       (grid[k]*2.0));
			   }
			}else{
			   if (i_ch->co[i][j]<0.0){
			      i_ch->co[i][j]=0.0001;
			   }else{
			      i_ch->co[i][j]*=2.0;
			   }
			}
		}
	}
 }
 for(i=1;i<=rows;i++){
	for(j=1;j<=cols;j++){
       if (Z0->co[i][j]==(double)UV->V->co[2]){
	      i_ch->co[i][j]=UV->V->co[2];
	   }
	}
 }

}

//***************************************************************************


//Distance from channel (in number of pixels)
void distance_from_channel(DOUBLEMATRIX *dist, SHORTMATRIX *DD, SHORTMATRIX *ST){

	static short r_DD[11]={0,0,-1,-1,-1,0,1,1,1,-9999,0}; /* differential of number-pixel for rows and*/
	static short c_DD[11]={0,1,1,0,-1,-1,-1,0,1,-9999,0}; /* columns, depending on Drainage Directions*/
	double b[10];                                         /* area perpendicular to the superficial flow divided by hsup*/
	long r, c, R, C;
	short dd;
 
	b[1]=0.;  
	b[2]=1.;             
	b[3]=sqrt(2.);  
	b[4]=1.;            
	b[5]=sqrt(2.);  
	b[6]=1.;  
	b[7]=sqrt(2.);   
	b[8]=1.;             
	b[9]=sqrt(2.); 
	
	initialize_doublematrix(dist, -9999.0);
	
	for(r=1;r<=dist->nrh;r++){
		for(c=1;c<=dist->nch;c++){
			if(DD->co[r][c]!=9){
																		
				dist->co[r][c]=0.0;
				
				if(ST->co[r][c]!=10){
					R=r;
					C=c;
					do{
						dd=DD->co[R][C];
						dist->co[r][c]+=b[dd+1];
						R+=r_DD[dd];
						C+=c_DD[dd];
					}while(ST->co[R][C]!=10 && DD->co[R][C]!=9);
				}
				
				if(DD->co[R][C]==9) dist->co[r][c]=-9999.0;
			}
						
		}
	}
}
	
//***************************************************************************
	
void distance_from_channel2(DOUBLEMATRIX *dist, SHORTMATRIX *ST, LONGVECTOR *rch, LONGVECTOR *cch){

	long r, c, R, C, ch;
	double d;

	initialize_doublematrix(dist, -9999.0);

	for(r=1;r<=dist->nrh;r++){
		for(c=1;c<=dist->nch;c++){
			if(ST->co[r][c]==10){
				dist->co[r][c]=0.0;
			}else if(ST->co[r][c]==0){
				dist->co[r][c]=1.E10;
				for(ch=1;ch<=rch->nh;ch++){
					R=rch->co[ch];
					C=cch->co[ch];
					d=sqrt(pow(R-r,2.)+pow(C-c,2.));
					if(d<dist->co[r][c]) dist->co[r][c]=d;
				}
			}
		}
	}
}

//***************************************************************************

void set_boundary_condition(DOUBLEMATRIX *Z, DOUBLEMATRIX *LC, short code, SHORTMATRIX *pixel_type, double novalue){
	
	//calculates the outlet of a basin, and assign value "code" to the matrix pixel_type to the outlet
	
	long r, c, rmin=0, cmin=0;
	long nr=Z->nrh, nc=Z->nch;
	double min=1.E99;
	
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			if(LC->co[r][c] != novalue){
				pixel_type->co[r][c]=0;
				if( Z->co[r][c] < min && (LC->co[r+1][c] != novalue || LC->co[r-1][c] != novalue || LC->co[r][c+1] != novalue
					|| LC->co[r][c-1] != novalue ) ){
					
					min = Z->co[r][c];
					rmin = r;
					cmin = c;
				}
			}
		}
	}
	
	if(rmin == 0 || cmin==0) t_error("Error in set_boundary_condition");
	
	pixel_type->co[rmin][cmin] = code;

}

//***************************************************************************

