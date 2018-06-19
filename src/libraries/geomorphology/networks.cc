#include "turtle.h"
#include "math.h"
#include "networks.h"
#include "t_utilities.h"
#include "t_datamanipulation.h"
#include "t_random.h"

#define CILINDRICAL 0

#define SQRT2  1.414213562373095
double weight[12]= {0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2,0,0,0};

short BOUNDARY=0;



/*----------------------------------------------------------------------------------*/

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




