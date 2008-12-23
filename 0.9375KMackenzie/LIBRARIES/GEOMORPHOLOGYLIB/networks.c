#include "turtle.h"
#include "math.h"
#include "networks.h"
#include "t_utilities.h"
#include "t_datamanipulation.h"
#include "t_random.h"

#define CILINDRICAL 0

#define SQRT2  1.414213562373095
double weight[12]={0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2,0,0,0};
 
short BOUNDARY=0;

/*--------------------------------------------------------------------------*/ 

LONGPAIR* initialize_flowdirections_with_outlets(SHORTMATRIX *m)

{



short u=1;

long i,j;

LONGPAIR *out=NULL,*tmp=NULL;

i=1;

for(j=1;j<=m->nch;j++){

	if(m->co[i][j]!=9){

		m->co[i][j]=9;

		if(u){

			printf("Warning::Network had no proper boundaries");

			u=0;

			

		}	

	}

}



i=m->nrh;



for(j=1;j<=m->nch;j++){

	if(m->co[i][j]!=9){

		m->co[i][j]=9;

		if(u){

			printf("Warning::Network had no proper boundaries");

			u=0;

			

		}	

	}

}





j=1;



for(i=2;i<=m->nrh-1;i++){

	if(m->co[i][j]!=9){

		m->co[i][j]=9;

		if(u){

			printf("Warning::Network had no proper boundaries");

			u=0;

			

		}	

	}

}



j=m->nch;



for(i=2;i<=m->nrh-1;i++){

	if(m->co[i][j]!=9){

		m->co[i][j]=9;

		if(u){

			printf("Warning::Network had no proper boundaries");

			u=0;

			

		}	

	}

}



for(i=2;i<=m->nrh-1;i++){

	for(j=2;j<=m->nch-1;j++){

		if(m->co[i][j]!=9){

			if(m->co[i][j]==10){

				tmp=new_longpair();

				out=appendto(out,tmp);

				tmp->i=i;

				tmp->j=j;

			}else{

				m->co[i][j]=0;

			}

		}

	}

}

	

return out;



}


void drainagedirections(DOUBLEMATRIX *elevations,SHORTMATRIX *directions,DOUBLEVECTOR *U,DOUBLEVECTOR *V)

{


long i,j,k;
static long     dir[11][3] =      {

																{ 0, 0, 0},

																{ 0, 1, 1},

																{-1, 1, 2},

																{-1, 0, 3},

																{-1,-1, 4},

																{ 0,-1, 5},

																{ 1,-1, 6},

																{ 1, 0, 7},

																{ 1, 1, 8},
																
																{ 0, 0, 9},

																{ 0, 0, 10}



													 };

double grid[9];

double df=0,steepestdescent=0;

/* char ch; */

grid[0]=1;
grid[1]=grid[5]=U->co[1];
grid[3]=grid[7]=U->co[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

initialize_shortmatrix(directions,9);

for(i=2;i<=elevations->nrh-1;i++){
	for(j=2;j<=elevations->nch-1; j++){
		if(elevations->co[i][j]!=V->co[2]){
			steepestdescent=0;
			for(k=1;k<=8;k++){
				if(elevations->co[i+dir[k][0]][j+dir[k][1]]!=V->co[2]){
						        
					df=(elevations->co[i][j]-elevations->co[i+dir[k][0]][j+dir[k][1]])/grid[k];
				  /*  printf("%d %d %d %d %f %f %f %f %d\n",i,j,i+dir[k][0],j+dir[k][1],elevations->co[i][j],elevations->co[i+dir[k][0]][j+dir[k][1]],df,steepestdescent,k);
				    scanf("%c",&ch); */

					if(df> steepestdescent){ 
						steepestdescent=df;
						directions->co[i][j]=k;
					}
				}
		    }
		    
			/*  printf("%d\n",directions->co[i][j]);  */
			if(directions->co[i][j]==9){
				if (is_ontheborder(elevations,V,i,j)){
					directions->co[i][j]=10;
			    }else {
					printf("Warning::a pit was found at position %d %d\n", i,j);
					directions->co[i][j]=12;
		
			    } 
				

		   }
	}
 }

}

}

/*--------------------------------------------------*/

void drainagedirections_modify(DOUBLEMATRIX *elevations,SHORTMATRIX *directions,DOUBLEVECTOR *U,DOUBLEVECTOR *V,
                                  SHORTMATRIX *segna)
/*
Inputs: 	1) Z0 matrice delle elevazioni
			2) segna matrice di controllo (segna dove è lago)
			3) U vettore che definisce la dimensione dei pixel
			4) V vettore che definisce i novalues
Outputs:	1) directions shortmatrix of directions
			2) segna matrice di controllo */
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
{
long i,j,k;
static long dir[11][3] = {
{ 0, 0, 0},
{ 0, 1, 1},
{-1, 1, 2},
{-1, 0, 3},
{-1,-1, 4},
{ 0,-1, 5},
{ 1,-1, 6},
{ 1, 0, 7},
{ 1, 1, 8},
{ 0, 0, 9},
{ 0, 0, 10}
};
double grid[9];
double df=0,steepestdescent=0;
short cont=0;
/*char ch;*/ 

grid[0]=1;
grid[1]=grid[5]=U->co[1];
grid[3]=grid[7]=U->co[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

initialize_shortmatrix(directions,9);
for(i=2;i<=elevations->nrh-1;i++){
	for(j=2;j<=elevations->nch-1; j++){
		if(elevations->co[i][j]!=V->co[2] && segna->co[i][j]!=1){
			steepestdescent=0;
			for(k=1;k<=8;k++){
				if(elevations->co[i+dir[k][0]][j+dir[k][1]]!=V->co[2]){	        
					df=(elevations->co[i][j]-elevations->co[i+dir[k][0]][j+dir[k][1]])/grid[k];
				    /*printf("%d %d %d %d %f %f %f %f %d\n",i,j,i+dir[k][0],j+dir[k][1],elevations->co[i][j],elevations->co[i+dir[k][0]][j+dir[k][1]],df,steepestdescent,k);
				    scanf("%c",&ch);*/
					if(df> steepestdescent){ 
						steepestdescent=df;
						directions->co[i][j]=k;
					}
				}
			}  
			if(directions->co[i][j]==9){
			/** chiama is_ontheborder(elevations,V,i,j) */ 
				if(is_ontheborder(elevations,V,i,j)){
					directions->co[i][j]=10;
					printf("Outlet at position %d %d\n", i,j);
					cont++;
			    }else{
					printf("Warning::a pit was found at position %d %d\n", i,j);
					directions->co[i][j]=12;
			    } 
		   	}
		}
		if(segna->co[i][j]==1)directions->co[i][j]=11;
 	}
}
if(cont>1) printf("Warning::this matrix has multiple outlets \n");
}

/*--------------------------------------------------*/

void drainagedirections_modify_f(DOUBLEMATRIX *elevations,SHORTMATRIX *directions,FLOATVECTOR *U,FLOATVECTOR *V,
                                  SHORTMATRIX *segna)
/*
Inputs: 	1) Z0 matrice delle elevazioni
			2) segna matrice di controllo (segna dove è lago)
			3) U vettore che definisce la dimensione dei pixel
			4) V vettore che definisce i novalues
Outputs:	1) directions shortmatrix of directions
			2) segna matrice di controllo */
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
{
long i,j,k;
static long dir[11][3] = {
{ 0, 0, 0},
{ 0, 1, 1},
{-1, 1, 2},
{-1, 0, 3},
{-1,-1, 4},
{ 0,-1, 5},
{ 1,-1, 6},
{ 1, 0, 7},
{ 1, 1, 8},
{ 0, 0, 9},
{ 0, 0, 10}
};
double grid[9];
double df=0,steepestdescent=0;
short cont=0;
/*char ch;*/ 

grid[0]=1;
grid[1]=grid[5]=U->co[1];
grid[3]=grid[7]=U->co[2];
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(grid[1]*grid[1]+grid[3]*grid[3]);

initialize_shortmatrix(directions,9);
for(i=2;i<=elevations->nrh-1;i++){
	for(j=2;j<=elevations->nch-1; j++){
		if(elevations->co[i][j]!=V->co[2] && segna->co[i][j]!=1){
			steepestdescent=0;
			for(k=1;k<=8;k++){
				if(elevations->co[i+dir[k][0]][j+dir[k][1]]!=V->co[2]){	        
					df=(elevations->co[i][j]-elevations->co[i+dir[k][0]][j+dir[k][1]])/grid[k];
				    /*printf("%d %d %d %d %f %f %f %f %d\n",i,j,i+dir[k][0],j+dir[k][1],elevations->co[i][j],elevations->co[i+dir[k][0]][j+dir[k][1]],df,steepestdescent,k);
				    scanf("%c",&ch);*/
					if(df> steepestdescent){ 
						steepestdescent=df;
						directions->co[i][j]=k;
					}
				}
			}  
			if(directions->co[i][j]==9){
			/** chiama is_ontheborder(elevations,V,i,j) */ 
				if(is_ontheborder_f(elevations,V,i,j)){
					directions->co[i][j]=10;
					printf("Outlet at position %d %d\n", i,j);
					cont++;
			    }else{
					printf("Warning::a pit was found at position %d %d\n", i,j);
					directions->co[i][j]=12;
			    } 
		   	}
		}
		if(segna->co[i][j]==1)directions->co[i][j]=11;
 	}
}
if(cont>1) printf("Warning::this matrix has multiple outlets \n");
}



/*--------------------------------------------------*/

short is_ontheborder(DOUBLEMATRIX *elevations,DOUBLEVECTOR *V,long i, long j)

/* serve per trovare i punti di uscita ed i pits */

{
short k;
static long dir[11][3] ={
{ 0, 0, 0},
{ 0, 1, 1},
{-1, 1, 2},
{-1, 0, 3},
{-1,-1, 4},
{ 0,-1, 5},
{ 1,-1, 6},
{ 1, 0, 7},
{ 1, 1, 8},
{ 0, 0, 9},
{ 0, 0, 10}
};

for(k=1;k<=8;k++){
	if(elevations->co[i+dir[k][0]][j+dir[k][1]]==V->co[2]){
		return 1;
	}
}
return 0;
}

/*--------------------------------------------------*/

short is_ontheborder_f(DOUBLEMATRIX *elevations,FLOATVECTOR *V,long i, long j)

/* serve per trovare i punti di uscita ed i pits */

{
short k;
static long dir[11][3] ={
{ 0, 0, 0},
{ 0, 1, 1},
{-1, 1, 2},
{-1, 0, 3},
{-1,-1, 4},
{ 0,-1, 5},
{ 1,-1, 6},
{ 1, 0, 7},
{ 1, 1, 8},
{ 0, 0, 9},
{ 0, 0, 10}
};

for(k=1;k<=8;k++){
	if(elevations->co[i+dir[k][0]][j+dir[k][1]]==V->co[2]){
		return 1;
	}
}
return 0;
}


/*---------------------------------------------------------------*/

void select_channels(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,DOUBLEMATRIX *curv, 
                     SHORTMATRIX *m,double threshold)

{


long i,j,flw[2],counter=0;
double laplacian,slope,area;

SHORTMATRIX *flow;

flow=new_shortmatrix(m->nrh,m->nch);
copy_shortmatrix(m,flow);

printf("threshold(%f)\n",threshold);
for(i=1;i<=m->nrh;i++){
	for(j=1;j<=m->nch;j++){
		if (m->co[i][j]!=9){
           laplacian=curv->co[i][j];
           slope=drainagedirections->co[i][j];
           area=ca->co[i][j];
		   if (sqrt(area)*slope < threshold || laplacian < 0){
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
			m->co[i][j]=flow->co[i][j];
			/**/
			go_downstream(flw,m->co[flw[0]][flw[1]],m->nch);		
			while( flow->co[flw[0]][flw[1]] < 9 && m->co[flw[0]][flw[1]]==0){
			   counter++;
			   m->co[flw[0]][flw[1]]=flow->co[flw[0]][flw[1]];
			   /*printf("%d %d %d\n", m->co[flw[0]][flw[1]],flw[0],flw[1]);*/
			   go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);
			   /* scanf("%c",&ch); */
				}
			}
		}
	}
free_shortmatrix(flow);
printf("Channel pixels are: %d\n",counter);
}


/*-----------------------------------------------------*/

void select_hillslopes(LONGMATRIX *ca,DOUBLEMATRIX *drainagedirections,DOUBLEMATRIX *curv, 
                       SHORTMATRIX *m,double threshold)

{
long i,j,flw[2],counter=0,hillslopecounter=0,area,max;
double laplacian,slope,valore;
SHORTMATRIX *flow;
short mode;

max=0;
mode=1;
printf("threshold(%f)\n",threshold);
flow=new_shortmatrix(m->nrh,m->nch);
copy_shortmatrix(m,flow);
	for(i=2;i<=m->nrh-1;i++){
		for(j=2;j<=m->nch-1;j++){
            		if(m->co[i][j]<9){
                                 laplacian=curv->co[i][j];
                                 slope=drainagedirections->co[i][j];
                                 area=ca->co[i][j];
                                 valore=sqrt((double)area)*(double)slope;
								if(valore < threshold || laplacian < 0){
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

printf("Channel pixels are: %d\n",counter);

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
			    	/*printf("%d %d %d\n", m->co[flw[0]][flw[1]],flw[0],flw[1]);*/
					go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);
					/* scanf("%c",&ch); */
				
				}
			
		
			}
			
		}
	}
free_shortmatrix(flow);
printf("Hillslope pixels are: %d\n",hillslopecounter);
}

/*-----------------------------------------------------*/

void select_hillslopes_f(LONGMATRIX *ca,FLOATMATRIX *drainagedirections,FLOATMATRIX *curv, 
                       SHORTMATRIX *m,double threshold)

{
long i,j,flw[2],counter=0,hillslopecounter=0,area,max;
double laplacian,slope,valore;
SHORTMATRIX *flow;
short mode;

max=0;
mode=1;
printf("threshold(%f)\n",threshold);
flow=new_shortmatrix(m->nrh,m->nch);
copy_shortmatrix(m,flow);
	for(i=2;i<=m->nrh-1;i++){
		for(j=2;j<=m->nch-1;j++){
            		if(m->co[i][j]<9){
                                 laplacian=curv->co[i][j];
                                 slope=drainagedirections->co[i][j];
                                 area=ca->co[i][j];
                                 valore=sqrt((double)area)*(double)slope;
								if(valore < threshold || laplacian < 0){
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

printf("Channel pixels are: %d\n",counter);

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
			    	/*printf("%d %d %d\n", m->co[flw[0]][flw[1]],flw[0],flw[1]);*/
					go_downstream(flw,flow->co[flw[0]][flw[1]],m->nch);
					/* scanf("%c",&ch); */
				
				}
			
		
			}
			
		}
	}
free_shortmatrix(flow);
printf("Hillslope pixels are: %d\n",hillslopecounter);
}



/*---------------------------------------------------------------------*/

void  outletdistance(SHORTMATRIX *m, DOUBLEMATRIX *dist, DOUBLEVECTOR *U)

{

long i,j,oldir,flow[2];
double grid[11],count,dx,dy;


dx=U->co[1];
dy=U->co[2];
 
grid[0]=grid[9]=grid[10]=0;
grid[1]=grid[5]=dx;
grid[3]=grid[7]=dy;
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(dx*dx+dy*dy);

 

for(i=1;i<=m->nrh;i++) {
	for(j=1;j<=m->nch;j++){
	        flow[0]=i;
		flow[1]=j;
		if(sourcesq(m,flow)==1){
				count=0;
/*				printf("#i=%d j=%d %f\n",i, j,count);			*/
 				oldir=m->co[flow[0]][flow[1]];
                go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
/*                printf("(%d,%d,%d,%f)\n",flow[0],flow[1],m->co[flow[0]][flow[1]],dist->co[flow[0]][flow[1]]);*/
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
				 	count+=grid[oldir];
/*		            printf("(%d,%d,%d,%f)\n",flow[0],flow[1],m->co[flow[0]][flow[1]],dist->co[flow[0]][flow[1]]);
				 	printf("%f %d \n",count,oldir); */
					oldir=m->co[flow[0]][flow[1]];	
                    go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(dist->co[flow[0]][flow[1]]> 0 ){
 				    count+=grid[oldir]+dist->co[flow[0]][flow[1]];
/* 				 	printf(" %d %d Sum: %f\n",i,j,count); */
 	                 			dist->co[i][j]=count;
 				}else if(m->co[flow[0]][flow[1]]>9){
 				       	dist->co[flow[0]][flow[1]]=0;
 					count+=grid[oldir];
/* 				 	printf(" %d %d Sum: %f\n",i,j,count); 	*/			
 					dist->co[i][j]=count;
 				}
/*	            scanf("%c",&ch); */
	            
 		        flow[0]=i;
 		        flow[1]=j;
                oldir=m->co[flow[0]][flow[1]];
 			/* count=dist->co[i][j]; */
 		        go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<0){
				 	count-=grid[oldir];
					dist->co[flow[0]][flow[1]]=count;
/* 				 	printf(" %d %d Sum: %f\n",flow[0],flow[1],dist->co[flow[0]][flow[1]]);
  				    scanf("%c",&ch); */
 					oldir=m->co[flow[0]][flow[1]];	
                    go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 		    
 		    
 		    }
 	   	
	}
/* printf("#%d/%d\t",i,m->nrh); */
}	

}


/*---------------------------------------------------------------------*/

void topological_outletdistance(SHORTMATRIX *m,LONGMATRIX *dist)

{

long i,j,oldir,flow[2];
double count;


 

for(i=1;i<=m->nrh;i++) {
	for(j=1;j<=m->nch;j++){
	        flow[0]=i;
		flow[1]=j;
		if(sourcesq(m,flow)==1){
				count=0;
/*				printf("#%d %d %f\n",i, j,count);			*/
 				oldir=m->co[flow[0]][flow[1]];
                go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
				 	count+=1;
/*				 	printf("%f %d \n",count,oldir); */
					oldir=m->co[flow[0]][flow[1]];	
                    go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(dist->co[flow[0]][flow[1]]> 0){
 				     	count+=1+dist->co[flow[0]][flow[1]];
/* 				 	printf(" %d %d Sum: %f\n",i,j,count); */
 	                 			dist->co[i][j]=count;
 				}else if(m->co[flow[0]][flow[1]]>9){
 				       	dist->co[flow[0]][flow[1]]=0;
 					count+=1;
/* 				 	printf(" %d %d Sum: %f\n",i,j,count); 	*/				
 					dist->co[i][j]=count;
 				}
/*	            scanf("%c",&ch); */
	            
 		        flow[0]=i;
 		        flow[1]=j;
                oldir=m->co[flow[0]][flow[1]];
 			/* count=dist->co[i][j]; */
 		        go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<0){
				 	count-=1;
					dist->co[flow[0]][flow[1]]=count;
 /* 				 	printf(" %d %d Sum: %f\n",flow[0],flow[1],dist->co[flow[0]][flow[1]]);
  				        scanf("%c",&ch); 
*/ 					oldir=m->co[flow[0]][flow[1]];	
                    go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 		    
 		    
 		    }
 	   	
	}
/* printf("#%d/%d\t",i,m->nrh); */
}	

}

/*-----------------------------------------------------------------------------------*/
void topological_hillslopes_channels_outletdistance(SHORTMATRIX *m,SHORTMATRIX* hillslopes, LONGMATRIX *dist)

{

long i,j,oldir,flow[2],cflow[2];
double counthillslope,count;

 

for(i=1;i<=m->nrh;i++) {
	for(j=1;j<=m->nch;j++){
	    flow[0]=i;
		flow[1]=j;
		if(sourcesq(m,flow)==1){
				counthillslope=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 &&
				                                      hillslopes->co[flow[0]][flow[1]] <9){
				 	counthillslope++;
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(hillslopes->co[flow[0]][flow[1]] >9){
 					counthillslope++;
 					dist->co[i][j]=-counthillslope;
 				} else if(dist->co[flow[0]][flow[1]] < 0){
 				    	counthillslope+=1-dist->co[flow[0]][flow[1]];
  					dist->co[i][j]=-counthillslope;
					while( m->co[flow[0]][flow[1]] <9 && 
				                                      hillslopes->co[flow[0]][flow[1]] <9){
                    					go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 					}
 				} else if(dist->co[flow[0]][flow[1]] == 0){
 				    	counthillslope++;
  					dist->co[i][j]=-counthillslope;
 				}	
 				
 				cflow[0]=flow[0];
 				cflow[1]=flow[1]; 				
	 		        	flow[0]=i;
 		        		flow[1]=j;
                			oldir=m->co[flow[0]][flow[1]];
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(  m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 && 
				                                       hillslopes->co[flow[0]][flow[1]] <9){
				 	counthillslope--;
					dist->co[flow[0]][flow[1]]=-counthillslope;
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
                
 		  /* */
 		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];
				count=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
				 	count++;
					oldir=m->co[flow[0]][flow[1]];		 			    
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(dist->co[flow[0]][flow[1]]>0){
 				     	count+=1+dist->co[flow[0]][flow[1]];
 	                 			dist->co[cflow[0]][cflow[1]]=count;
 				}else if(m->co[flow[0]][flow[1]] > 9){
  					count++;
 					dist->co[cflow[0]][cflow[1]]=count;
 				}
	            

		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];
                			oldir=m->co[flow[0]][flow[1]];
 				/* count=dist->co[i][j]; */
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(   m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
 					count--;
 					dist->co[flow[0]][flow[1]]=count;
  		      		  	go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
					  	
 					
 				}
 		 /* */     
 		    }
 	   	
	}
  
}

}

/*-----------------------------------------------------*/

void hillslopes_channels_outletdistance(SHORTMATRIX *m,SHORTMATRIX* hillslopes, DOUBLEMATRIX *dist,DOUBLEVECTOR *U)

{


long i,j,oldir,flow[2],cflow[2],max;
double grid[11],counthillslope,count,dx,dy;
max=0;
dx=U->co[1];
dy=U->co[2];
grid[0]=grid[9]=grid[10]=0;
grid[1]=grid[5]=dx;
grid[3]=grid[7]=dy;
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(dx*dx+dy*dy);
for(i=2;i<=m->nrh-1;i++) {
	for(j=2;j<=m->nch-1;j++){
	    flow[0]=i;
		flow[1]=j;
		if(sourcesq(m,flow)==1){
				counthillslope=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 &&
				                                      hillslopes->co[flow[0]][flow[1]] <9){
				 	counthillslope+=grid[oldir];
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(hillslopes->co[flow[0]][flow[1]] >9){
 					counthillslope+=grid[oldir];
 					dist->co[i][j]=-counthillslope;
 				} else if(dist->co[flow[0]][flow[1]] < 0){
 				    	counthillslope+=grid[oldir]-dist->co[flow[0]][flow[1]];
  					dist->co[i][j]=-counthillslope;
					while( m->co[flow[0]][flow[1]] <9 && 
				          		hillslopes->co[flow[0]][flow[1]] <9){
                    					go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 					}
					
 				} else if(dist->co[flow[0]][flow[1]] == 0){
 				    	counthillslope+=grid[oldir];
  					dist->co[i][j]=-counthillslope;
 				}	
 				
 				cflow[0]=flow[0];
 				cflow[1]=flow[1];
    
     
		 		flow[0]=i;
 		        		flow[1]=j;
                			oldir=m->co[flow[0]][flow[1]];
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(  m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 && 
				                                       hillslopes->co[flow[0]][flow[1]] <9){
					counthillslope-=grid[oldir];
					dist->co[flow[0]][flow[1]]=-counthillslope;
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
                
 		  /* */

 		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];

				count=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
				 	count+=grid[oldir];
					oldir=m->co[flow[0]][flow[1]];	 			    
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(dist->co[flow[0]][flow[1]]>0){
 					count+=grid[oldir]+dist->co[flow[0]][flow[1]];
 	                 			dist->co[cflow[0]][cflow[1]]=count;
 				}else if(m->co[flow[0]][flow[1]] > 9){
  					count+=grid[oldir];
 					dist->co[cflow[0]][cflow[1]]=count;
 				}
	            

		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];
                			oldir=m->co[flow[0]][flow[1]];
 				/* count=dist->co[i][j]; */
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
 					count-=grid[oldir];
 					oldir=m->co[flow[0]][flow[1]];
 					dist->co[flow[0]][flow[1]]=count;
  		      		  	go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				}
 		 /* */     
 		    }
 	   	
	}
  
}

}

/*-----------------------------------------------------*/

void hillslopes_channels_outletdistance_f(SHORTMATRIX *m,SHORTMATRIX* hillslopes, DOUBLEMATRIX *dist,FLOATVECTOR *U)

{


long i,j,oldir,flow[2],cflow[2],max;
double grid[11],counthillslope,count,dx,dy;
max=0;
dx=U->co[1];
dy=U->co[2];
grid[0]=grid[9]=grid[10]=0;
grid[1]=grid[5]=dx;
grid[3]=grid[7]=dy;
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(dx*dx+dy*dy);
for(i=2;i<=m->nrh-1;i++) {
	for(j=2;j<=m->nch-1;j++){
	    flow[0]=i;
		flow[1]=j;
		if(sourcesq(m,flow)==1){
				counthillslope=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 &&
				                                      hillslopes->co[flow[0]][flow[1]] <9){
				 	counthillslope+=grid[oldir];
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(hillslopes->co[flow[0]][flow[1]] >9){
 					counthillslope+=grid[oldir];
 					dist->co[i][j]=-counthillslope;
 				} else if(dist->co[flow[0]][flow[1]] < 0){
 				    	counthillslope+=grid[oldir]-dist->co[flow[0]][flow[1]];
  					dist->co[i][j]=-counthillslope;
					while( m->co[flow[0]][flow[1]] <9 && 
				          		hillslopes->co[flow[0]][flow[1]] <9){
                    					go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 					}
					
 				} else if(dist->co[flow[0]][flow[1]] == 0){
 				    	counthillslope+=grid[oldir];
  					dist->co[i][j]=-counthillslope;
 				}	
 				
 				cflow[0]=flow[0];
 				cflow[1]=flow[1];
    
     
		 		flow[0]=i;
 		        		flow[1]=j;
                			oldir=m->co[flow[0]][flow[1]];
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(  m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]==0 && 
				                                       hillslopes->co[flow[0]][flow[1]] <9){
					counthillslope-=grid[oldir];
					dist->co[flow[0]][flow[1]]=-counthillslope;
					oldir=m->co[flow[0]][flow[1]];	
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
                
 		  /* */

 		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];

				count=0;
 				oldir=m->co[flow[0]][flow[1]];
                			go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while( m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
				 	count+=grid[oldir];
					oldir=m->co[flow[0]][flow[1]];	 			    
                    				go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
 				}
 				if(dist->co[flow[0]][flow[1]]>0){
 					count+=grid[oldir]+dist->co[flow[0]][flow[1]];
 	                 			dist->co[cflow[0]][cflow[1]]=count;
 				}else if(m->co[flow[0]][flow[1]] > 9){
  					count+=grid[oldir];
 					dist->co[cflow[0]][cflow[1]]=count;
 				}
	            

		        		flow[0]=cflow[0];
 		        		flow[1]=cflow[1];
                			oldir=m->co[flow[0]][flow[1]];
 				/* count=dist->co[i][j]; */
 		        		go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				while(m->co[flow[0]][flow[1]] <9 && dist->co[flow[0]][flow[1]]<=0){
 					count-=grid[oldir];
 					oldir=m->co[flow[0]][flow[1]];
 					dist->co[flow[0]][flow[1]]=count;
  		      		  	go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
				}
 		 /* */     
 		    }
 	   	
	}
  
}

}



/*----------------------------------------------------------------------------------*/ 
/* 	Inputs:  1) m shortmatrix of flowing directions
	Outputs: 1) ca matrix that will contain the values of contributing area [pixels]*/
void tca(SHORTMATRIX *flow,LONGMATRIX *ca)
{
/* char ch;*/
long i,j,punto[2];
//long yop=0;

if(flow->nrh !=ca->nrh || flow->nch !=ca->nch){
	t_error("Flow matrix and contributing area matrix have not the same dimensions");
}
for(i=1;i<=ca->nrh;i++) {
	meter(i,ca->nrh,10,"Lines of matrix ca parsed","\n");
	for(j=1;j<=ca->nch;j++) {
		if(flow->co[i][j]!=9 ){
	    	punto[0]=i;punto[1]=j;
	    	while(flow->co[punto[0]][punto[1]]<9){
				(ca->co[punto[0]][punto[1]])++;
				
				/* lancia go_downstream 
					Inputs: 1) the position; 
							2) the flowing direction (coded); 
							3) the number of columns in matrix 
							(It is not used if the constant CILINDRICAL is not set to 1) */
				
				go_downstream(punto,flow->co[punto[0]][punto[1]],ca->nch);
				
			}
			if(flow->co[punto[0]][punto[1]]==10)(ca->co[punto[0]][punto[1]])++;
		}		
	}
}
}

/*-----------------------------------------------------------------------------*/





/*-----------------------------------------------------------------------------*/
/* commented because the Warning: variable 'work' is not initialized before being used
void randomflow(long x,long y, SHORTMATRIX *m)
{
long r,k=0,p;
LONGPAIR *work,*home;
LONGPAIR *now;

if(m->co[x][y]!=11) t_error("Not a border point");

now=work;

r=m->co[x+work->i][y+work->j];

while(r==9 || r==0 || r==11 ||

      crossQ(x,y,x+now->i,y+now->j,m) ){

       if(work->next==NULL) {

            work=prependto(home,work);

       		home=NULL;

       		printf("\nWarning::pit in position {%d,%d}\n",x,y);

       		m->co[x][y]=12;

       		}

       else{

       		work=work->next;

       		now->next=NULL;

       		home=prependto(home,now);

	   		now=work;

       		r=m->co[x+work->i][y+work->j];
       		

       }

}



p=coords2int(x,y,x+now->i,y+now->j);

m->co[x][y]=p;

work=rotateleft(work);

work=appendto(work,home);

home=NULL;

 

}

*/



/*------------------------This is net.c------------------------*/







 

/* This make a point to be an outlet: 

the eight direction topology is assumed. */



/*

void outlet(unsigned int i,unsigned int j, umatrix *m)



{

if(i==0 || j==0 || i==((m->row)-1) || j==((m->col)-1)) 

ocnerror("The outlet cannot be allocated at this position\n");

else if( element(m,i-1,j-1)==0.) element(m,i,j)=10;

else if( element(m,i-1,j)==0.) element(m,i,j)=10;

else if( element(m,i-1,j+1)==0.) element(m,i,j)=10;

else if( element(m,i,j-1)==0.) element(m,i,j)=10;

else if( element(m,i,j)==0.) element(m,i,j)=10;

else if( element(m,i,j+1)==0.) element(m,i,j)=10;

else if( element(m,i+1,j-1)==0.) element(m,i,j)=10;

else if( element(m,i+1,j)==0.) element(m,i,j)=10;

else if( element(m,i+1,j+1)==0.) element(m,i,j)=10;

}



*/



/*--------------------------------------------------------------------------*/ 

/* Checking the correctness of the new link */



long crossQ(long  i,long j,long flowi,long flowj,SHORTMATRIX *m)

{

if((flowi-i)==-1 && (flowj-j)==1){

	if(m->co[i-1][j]==8 || m->co[i][j+1]==4) return 1;

	else return 0;  

	}

else if((flowi-i)==1 && (flowj-j)==1){

	if(m->co[i][j+1]==6 || m->co[i+1][j]==2) return 1;

	else return 0; 

	}

	

else if((flowi-i)==+1 && (flowj-j)==-1){

	if(m->co[i][j-1]==8 || m->co[i+1][j]==4) return 1;

	else return 0; 

	}

else if((flowi-i)==-1 && (flowj-j)==-1){

	if(m->co[i][j-1]==2 || m->co[i-1][j]==6) return 1;

	else return 0;  

	}

else return 0;



}





/*--------------------------------------------------------------------------*/ 

/* Transforming the coordinate differences between the point

and the  drained point into integers */



short coords2int(long i,long j,long flowi,long flowj)



{



short dify=flowi-i,difx=flowj-j;



if(difx==-1 && dify==-1)      return 4;

else if(difx==-1 && dify==0)  return 5;

else if(difx==-1 && dify==+1) return 6;

else if(difx== 0 && dify==-1) return 3;

else if(difx== 0 && dify==+1) return 7;

else if(difx==+1 && dify==-1) return 2;

else if(difx==+1 && dify== 0) return 1;

else if(difx==+1 && dify==+1) return 8;

else if(difx== 0 && dify== 0 && i!=-1 && j!=-1) return 0;

else  return 12;





}





/*--------------------------------------------------------------------------

go_downstream returns the row and column indexes corresponding to

the flow directions.

----------------------------------------------------------------------------*/


void go_downstream_2(long *punto,short number,long max)

{

 

 

 	   static int v[13][2] = {       

                             { 0, 0},

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

                             { 0, 0}
                         

                            };



number=abs(number);

 if (number < 0 || number > 12){

    t_error("Wrong number in flow matrix");

}



switch (number) {

                        case 11:

                                punto[0] = punto[1] = -1;

                                break;

                        

                        case 12:

                                punto[0] = punto[1] = -1;

                                break;
                                                        

						default:

                                punto[0] +=  v[number][0];

                                punto[1] +=  v[number][1];

                                break;

                }

       if(punto[1]==1){
       		
       		punto[1]=max-1;
       		BOUNDARY=1;
       		
       }else if (punto[1]==max){
            
            punto[1]=2;
       		BOUNDARY=1;

       }

}


void go_downstream(long *punto,short number,long max)

/* Inputs: 1) the position; 
		2) the flowing direction (coded); 
		3) the number of columns in matrix (It is not used if the constant CILINDRICAL is not set to 1) */

{
 	   static int v[13][2] = {       

                             { 0, 0},

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

                             { 0, 0}              
                            };

number=abs(number);
if (number < 0 || number > 12){
    t_error("Wrong number in flow matrix");
}

max=0;
switch (number) {
	case 12:
	punto[0] = punto[1] = -1;
	break;                                                 
	default:
	punto[0] +=  v[number][0];
	punto[1] +=  v[number][1];
	break;
}
}

/*
tcamax restituisce un numero: 0 nel caso che il pixel che presenta distanza dalla sorgente diss e area cumulata
maz rappresenta la direzione principale in una diramazione. 1 se rappresenta un affluente alla direzione principale
*/

long tcamax(SHORTMATRIX *m,LONGMATRIX *tca,DOUBLEMATRIX *dist,long *flow,long maz,double diss)

{

 static  short     dir[11][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0},
                                                { 0, 0, 0}
                                             
                                        };

long k;

for(k=1; k<=8; k++){
	if(m->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]==dir[k][2]) {
		if(tca->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]>=maz) {
			if(tca->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]==maz){ 
				if(dist->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]>diss) return  1;
			}
		else return 1;
		}
	}
}
return 0;


}  


/*--------------------------------------------------------------------------
NB. INSTEAD OF omatrix one could alocate the values in a
XY structure saving memory usage
----------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------
select_EOLink definisce la presenza di Link rilascaindo una matrice con tutti zero
a parte nei pixel che rappresentano la fine dei Link. Il valore e' il valore
della matrice imatrix che e' la matrice di input
Il programma necessita della matrice delle direzioni di drenaggio e della magnitudo.
Un Link finisce quando si ha la confluenza di due percorsi in sostanza quando 
il valore di magnitudo aumenta durante un percorso
--------------------------------------------------------------------------*/
/*long select_EOLink(double th,SHORTMATRIX *m,LONGMATRIX *magn,DOUBLEMATRIX *imatrix, 
                              DOUBLEMATRIX *omatrix)
{
long i,j,flow[2];
double valore, valore_prec;
long counter=0;
for(i=1;i<=m->nrh;i++){
for(j=1;j<=m->nch;j++){
if(m->co[i][j]<9 && magn->co[i][j]>th){
flow[0]=i;
flow[1]=j;
valore_prec=magn->co[flow[0]][flow[1]];
/*viene memorizzato il valore di magn. del pixel in questione*/
/*go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
/*viene memorizzato il valore di magn. del pixel di drenaggio*/
/*valore=magn->co[flow[0]][flow[1]];
if(valore_prec!=valore){
   omatrix->co[i][j]=imatrix->co[i][j];
   counter++;
}
}
if(m->co[i][j]==10 && magn->co[i][j]>th){
                  omatrix->co[i][j]=imatrix->co[i][j];
                  counter++;
}
}
}
return counter;
}
*/
/*--------------------------------------------------------------------------
go_upstream_a definisce il numero di pixel che drenano nel pixel
definito da *p. LA routine definisce con p di uscita il pixel di drenaggio 
che rappresenta il ramo princiapale. Nel caso vi siano piu' pixel drenanti
quello che definisce il ramo princiaple viene definito usando l'area cumulata e la lunghezza
dando prevalenza alla prima
----------------------------------------------------------------------------*/
long go_upstream_a(double th,long *p, long *kk,SHORTMATRIX *m,LONGMATRIX *tca,DOUBLEMATRIX *l)

{

long area=0, lenght=0,point[2]={0,0},count=0;

static  short     dir[13][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0},
                                                { 0, 0, 0},
                                                { 0, 0, 0},
                                                { 0, 0, 0} 
                                        };

register short k;
          point[0]=p[0];
          point[1]=p[1];
	for(k=1; k<=8; k++){
		if(abs(m->co[p[0]+dir[k][0]][p[1]+dir[k][1]])==dir[k][2]) {
			count++;
			if(tca->co[p[0]+dir[k][0]][p[1]+dir[k][1]]>=area && 
			   tca->co[p[0]+dir[k][0]][p[1]+dir[k][1]]>=th){
			           if(tca->co[p[0]+dir[k][0]][p[1]+dir[k][1]]==area){
			              if(l->co[p[0]+dir[k][0]][p[1]+dir[k][1]]>lenght){
                				*kk=k;
                				area=tca->co[p[0]+dir[k][0]][p[1]+dir[k][1]];
		       			lenght=l->co[p[0]+dir[k][0]][p[1]+dir[k][1]];
					point[0]=p[0]+dir[k][0];
					point[1]=p[1]+dir[k][1];
			              }
			           }else{
			               *kk=k;
			          	    area=tca->co[p[0]+dir[k][0]][p[1]+dir[k][1]];
				    lenght=l->co[p[0]+dir[k][0]][p[1]+dir[k][1]];
				    point[0]=p[0]+dir[k][0];
				    point[1]=p[1]+dir[k][1];
				}
			}
		}
	
	}
	p[0]=point[0];
	p[1]=point[1];
	
	return count;
}




/*--------------------------------------------------------------------------
UNDER COSTRUCTION
----------------------------------------------------------------------------
short go_upstream_b(long *p,SHORTMATRIX *m,LONGMATRIX *tca,DOUBLEMATRIX *length)

{

long area=0;
short point[2];

static  short     dir[10][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0} 
                                             
                                        };

register short k,count=0;

    
	for(k=1; k<=8; k++){
		if(abs(m->co[i+dir[k][0]][j+dir[k][1]])==dir[k][2]) {
			count++;
			if(tca->elment[i+dir[k][0]][j+dir[k][1]]>area){
				area=tca->co[i+dir[k][0]][j+dir[k][1]]);
				
				point[0]=i+dir[k][0];
				point[1]=j+dir[k][1];
			} else if(tca->elment[i+dir[k][0]][j+dir[k][1]]==area){
				if(element->length[i+dir[k][0]][j+dir[k][1]]>ll){
					
				}
			}
		}
	
	}
	
	p[0]=point[0];
	p[1]=point[1];
	return count;
}  


*/
/*--------------------------------------------------------------------------*/ 

/* Print the umatrix in Mathematica format in the specified file  */



/*--------------------------------------------------------------------------*/ 





double energyexpenditure(LONGMATRIX *tca,double ex,long th)



{



long i,j;

double enex=0;





if(th>=0){



	for(i=1;i<=tca->nrh;i++){

		for(j=2;j<=tca->nch-1;j++){

				if(tca->co[i][j] >th) {

					enex+=pow(tca->co[i][j],ex);

				}

			

		}

	/* printf("E=%f\n",enex); */



	}



} else {



	t_error("A negative threshold is not allowed");

}



return enex;



}


/*--------------------------------------------------------------------------*/ 





double r_energyexpenditure(DOUBLEMATRIX *tca,double ex,double th)



{



long i,j;

double enex=0;





if(th>=0){



	for(i=1;i<=tca->nrh;i++){

		for(j=2;j<=tca->nch-1;j++){

				if(tca->co[i][j] >th) {

					enex+=pow(tca->co[i][j],ex);

				}

			

		}

	/* printf("E=%f\n",enex); */



	}



} else {



	t_error("A negative threshold is not allowed");

}



return enex;



}


/*--------------------------------------------------------------------------*/ 





double weighted_energyexpenditure(SHORTMATRIX *m,LONGMATRIX *tca,double ex,long th)



{



long i,j;

double enex=0;

extern double weight[12];




if(th>=0){



	for(i=1;i<=tca->nrh;i++){

		for(j=1;j<=tca->nch;j++){

				if(tca->co[i][j] >th) {

					enex+=pow(tca->co[i][j],ex)*weight[m->co[i][j]];

				}

			

		}

	/* printf("E=%f\n",enex); */



	}



} else {



	t_error("A negative threshold is not allowed");

}



return enex;



}



/*--------------------------------------------------------------------------*/ 





double r_weighted_energyexpenditure(SHORTMATRIX *m,DOUBLEMATRIX *tca,double ex,double th)



{



long i,j;

double enex=0;

extern double weight[12];




if(th>=0){



	for(i=1;i<=tca->nrh;i++){

		for(j=1;j<=tca->nch;j++){

				if(tca->co[i][j] >th) {

					enex+=pow(tca->co[i][j],ex)*weight[m->co[i][j]];

				}

			

		}

	/* printf("E=%f\n",enex); */



	}



} else {



	t_error("A negative threshold is not allowed");

}



return enex;



}




/*--------------------------------------------------------------------------*/ 



LONGPAIR * topology(void )

{


LONGPAIR  NW,  N,  NE; 
LONGPAIR   W,  H,   E;
LONGPAIR  SW,  S,  SE;

long jj=1747,j;

LONGPAIR * work;







SW.i=W.i=NW.i=-1; 

N.i=S.i= H.i=0;  

NE.i=E.i=SE.i=1; 



SW.j=S.j=SE.j=+1; 

W.j=H.j=E.j=0; 

NW.j=N.j=NE.j=-1; 



work=&SW;





SW.next=&E;

E.next=&SE;

SE.next=&S;

S.next=&W;

W.next=&N;

N.next=&NW;

NW.next=&NE;

NE.next=NULL;  



j=floor(8*ran1(&jj));
work=rotate(work,j);



return  work;



}








/*--------------------------------------------------------------------------*/



void initialize_flowdirections(SHORTMATRIX* ca)



{



long i,j;



for(i=1;i<=ca->nrh;i++) {

	for(j=1;j<=ca->nch;j++){

		if(!(ca->co[i][j]==9 || ca->co[i][j]==10)){

			ca->co[i][j]=0; 

		}

	}

}	

}

/*-----------------------------------------------------------------------------------
sourcesq
----------------------------------------------------------s--------------------------*/

long sourcesq(SHORTMATRIX *m,long *flow)

{

 static  short     dir[11][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0},
                                                { 0, 0, 0}
                                             
                                        };

long k;

if(m->co[flow[0]][flow[1]] < 9 && m->co[flow[0]][flow[1]] > 0) {

	for(k=1; k<=8; k++){
		if(m->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]==dir[k][2]) {
			return 0;		
		}
	}
	return 1;
} else {

	return 0;
}

}

/*-----------------------------------------------------------------------------------
linkq
-------------------------------------------------------------------------------------*/

long linkq(SHORTMATRIX *m,long *flow)

{

 static  short     dir[11][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0},
                                                { 0, 0, 0}
                                             
                                        };

short k,counter=0;

if(m->co[flow[0]][flow[1]] <9 && m->co[flow[0]][flow[1]] >0) {

	for(k=1; k<=8; k++){
		if(m->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]==dir[k][2]) {
			counter++;
			if(counter >1){
				return 1;
			}		
		}
	}
	return 0;
} else {

	return 0;
}

}


/*---------------------------------------------------------------------*/

long conditioned_linkq(SHORTMATRIX *m,DOUBLEMATRIX *ca,long *flow,double th)

{

 static  short     dir[11][3] =      {       
                                                { 0, 0, 0},
                                                { 0, 1, 5},
                                                {-1, 1, 6},
                                                {-1, 0, 7},
                                                {-1,-1, 8},
                                                { 0,-1, 1},
                                                { 1,-1, 2},
                                                { 1, 0, 3},
                                                { 1, 1, 4},
                                                { 0, 0, 0},
                                                { 0, 0, 0}
                                             
                                        };

short k,counter=0;

if(m->co[flow[0]][flow[1]] <9 && m->co[flow[0]][flow[1]] >0) {

	for(k=1; k<=8; k++){
		if(m->co[flow[0]+dir[k][0]][flow[1]+dir[k][1]]==dir[k][2] && ca->co[flow[0]][flow[1]] >th) {
			counter++;
			if(counter >1){
				return 1;
			}		
		}
	}
	return 0;
} else {

	return 0;
}

}


/*---------------------------------------------------------------------*/

void downstream(SHORTMATRIX *m,DOUBLEMATRIX *quantity,DOUBLEMATRIX *lquantity,long lag,double dx, double dy,double novalue)

{

long i,j,oldir,oldirl,flow[2],flowl[2],oldflow[2];
double grid[11],count,reallag,lreallag;

if(lag <0){
	t_error("Lag value not allowed");
} else if(lag==0){
	copy_doublematrix(quantity,lquantity);
} else {
if(dx <dy ){
	reallag=dx*(lag+0.5);
	lreallag=reallag-dx;
}else{
	reallag=dy*(lag+0.5);
	lreallag=reallag-dy;

} 
grid[0]=grid[9]=grid[10]=0;
grid[1]=grid[5]=dx;
grid[3]=grid[7]=dy;
grid[2]=grid[4]=grid[6]=grid[8]=sqrt(dx*dx+dy*dy);
for(i=1;i<=quantity->nrh;i++) {
	for(j=1;j<=quantity->nch;j++){
	    flow[0]=i;
		flow[1]=j;
		flowl[0]=i;
		flowl[1]=j;
		
		if(sourcesq(m,flow)==1){
				count=0;
 				oldirl=m->co[flowl[0]][flowl[1]];
 				oldir=oldirl;
 				oldflow[0]=flowl[0];
 				oldflow[1]=flowl[1];
                go_downstream(flowl,m->co[flowl[0]][flowl[1]],m->nch);
                count+=grid[oldirl];
                while(count <= reallag && m->co[flowl[0]][flowl[1]]< 9 ){
                	oldirl=m->co[flowl[0]][flowl[1]];
        			oldflow[0]=flowl[0];
        			oldflow[1]=flowl[1];
                	go_downstream(flowl,m->co[flowl[0]][flowl[1]],m->nch);
           	        count+=grid[oldirl];

                }
                

                if(count-grid[oldirl]  > lreallag){
 				    lquantity->co[flow[0]][flow[1]]=quantity->co[oldflow[0]][oldflow[1]];
				} else {
 						lquantity->co[flow[0]][flow[1]]=novalue;
                }
 					
 				while(m->co[flowl[0]][flowl[1]] <9 ){
 					oldir=m->co[flow[0]][flow[1]];
 					count-=grid[oldir];
 					go_downstream(flow,m->co[flow[0]][flow[1]],m->nch);
					if(lquantity->co[flow[0]][flow[1]]!=novalue){ 
						break;
						}
 					while(count <= reallag && m->co[flowl[0]][flowl[1]]< 9 ){
                		oldirl=m->co[flowl[0]][flowl[1]];
        				oldflow[0]=flowl[0];
        				oldflow[1]=flowl[1];
                		go_downstream(flowl,m->co[flowl[0]][flowl[1]],m->nch);
           	        	count+=grid[oldirl];
                	}
	                if(count-grid[oldirl]  > lreallag){
 					    lquantity->co[flow[0]][flow[1]]=quantity->co[oldflow[0]][oldflow[1]];
        	        } else {
 							lquantity->co[flow[0]][flow[1]]=novalue;
                	}                	
 				
 				}
 				
 		}	
	            
 		    
 		
	}
}	
}
}

/*-------------------------------------------------------------------------*/

double network_dowstreamcorrelation(DOUBLEMATRIX* dist,DOUBLEMATRIX *ca,double mean1,double mean2,double novalue1,double novalue2)

{

double cr=0,dr=0,er=0,n=0;
long i,j;

if(dist==NULL || ca==NULL || dist->co==NULL || ca->co==NULL){
	t_error("Matrixes were not allocated");
} else if(dist->isdynamic!=1 || dist->isdynamic!=1 || dist->nrh <2 ){
    t_error("Matrixes were not properly allocated");
} else if(dist->nrh !=ca->nrh || dist->nch!=ca->nch) {
    t_error("These matrixes do not have the same dimensions");
}

for(i=1;i<=dist->nrh;i++){
	for(j=1;j<=dist->nch;j++){		
		if(dist->co[i][j]!=novalue1 && ca->co[i][j]!=novalue2){
			cr+=dist->co[i][j]*ca->co[i][j];
			dr+=dist->co[i][j]*dist->co[i][j];
			er+=ca->co[i][j]*ca->co[i][j];
			n++; 	
		}
	}
}

return (cr/n -mean1*mean2)/sqrt((dr/n-mean1*mean1)*(er/n-mean2*mean2));

}

/*---------------------------------------------------------------------*/

void sum_downstream(SHORTMATRIX *flow,DOUBLEMATRIX *la,DOUBLEMATRIX *dist)



{



/* char ch */

long i,j, punto[2];
long yop=0;



if(flow->nrh !=la->nrh || flow->nch !=la->nch || flow->nrh !=dist->nrh  || flow->nch !=dist->nch ){

	t_error("Flow matrix and the other matrixes do not have  the same dimensions");

}

for(i=1;i<=la->nrh;i++) {
	for(j=1;j<=la->nch;j++) {
		if(flow->co[i][j]!=9) {
			dist->co[i][j]=0;
}}}


for(i=1;i<=la->nrh;i++) {

    meter(i,la->nrh,10,"Lines of  matrix parsed","\n");	

	for(j=1;j<=la->nch;j++) {

		if(flow->co[i][j]!=9) {

	    	    punto[0]=i;punto[1]=j;

	    		while(flow->co[punto[0]][punto[1]] < 9){

	    		    /* printf("%d&",flow->co[punto[0]][punto[1]]);*/

					 dist->co[punto[0]][punto[1]]+=la->co[i][j];

					/* 

					printf("%d %d\n",punto[0],punto[1]);

					scanf("%c",&ch);

					*/

					go_downstream(punto,flow->co[punto[0]][punto[1]],la->nrh);

			}

	    dist->co[punto[0]][punto[1]]+=la->co[i][j];

		 

		}

	}

}



}


/*--------------------------------------------------------------------------------------*/

void dem_array_check(DOUBLEVECTOR *T, DOUBLEVECTOR *U )

{

 if(T->nh < 4 || U->nh <4  ) t_error("DEM header incorrect");
	if(T->co[1]!=U->co[1] || T->co[2]!=U->co[2] || T->co[4]!=U->co[4] || T->co[3]!=U->co[3]){
		printf("Header Elements:\n1:%f %f\n",T->co[1],U->co[1]);
		printf("2:%f %f\n",T->co[2],U->co[2]);
		printf("3:%f %f\n",T->co[3],U->co[3]);
		printf("4:%f %f\n",T->co[4],U->co[4]);
		
		t_error("DEM headers do not refer the same DEM");

  }
  
  
}
