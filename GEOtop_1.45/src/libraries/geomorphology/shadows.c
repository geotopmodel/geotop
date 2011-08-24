/*----------------------------------------------------------------------------------------------------------*/
/*Used to calculate shadow matrix*/
/*----------------------------------------------------------------------------------------------------------*/


#include "../fluidturtle/turtle.h"
#include "shadows.h"



void Orizzonte1(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{

/*=====================*/
long i,j,jj,y,I,J;
double zenith;
/*======================*/

if(beta!=0){
  	for(j=1;j<=quadrata;j++){
      	I=0;
       	J=0;
       	y=0;
       	for(jj=j;jj>0;jj--){
          	for(i=floor(1/tan(beta)*(j-jj))+1;
              	i<=floor(1/tan(beta)*(j-jj+1)) && i<=Z0->nrh;i++){
               	if(jj<=Z0->nch && Z0->co[i][jj]!=novalue){
				  	/*shadow->co[i][jj]=j;}}}}}*/
                  	if(curv->co[i][jj]==1 && I==0){
                     	I=i;
                     	J=jj;
                     	y=1;
                  	}else if(curv->co[i][jj]==1 && I!=0){
                   		zenith=(Z0->co[I][J]-Z0->co[i][jj])
                            /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-jj)*(double)delta,(double)2));
                    	if(zenith<=tan(alfa)){
                           	shadow->co[i][jj]=0;
                            I=i;
                            J=jj;
                        }else{
                             shadow->co[i][jj]=1;
                        }
                  	}else if(curv->co[i][jj]==0 && y==1){
                        zenith=(Z0->co[I][J]-Z0->co[i][jj])
                        	/sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-jj)*(double)delta,(double)2));
                   		if(zenith<=tan(alfa)){
                           	shadow->co[i][jj]=0;
                            y=0;
                        }else{
                             shadow->co[i][jj]=1;
                        }
                  	}
              	}
        	}
     	}
  	}
}else{
   	for(j=1;j<=Z0->nch;j++){
       	I=0;
       	J=0;
       	y=0;
       	for(i=1;i<=Z0->nrh;i++){
           	if(Z0->co[i][j]!=novalue){
              	if(curv->co[i][j]==1 && I==0){
                 	I=i;
                 	J=j;
                 	y=1;
              	}else if(curv->co[i][j]==1 && I!=0){
                    zenith=(Z0->co[I][J]-Z0->co[i][j])
                      	/sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                    if(zenith<=tan(alfa)){
                      	shadow->co[i][j]=0;
                       	I=i;
                       	J=j;
                    }else{
                       	shadow->co[i][j]=1;
                    }
              	}else if(curv->co[i][j]==0 && y==1){
                   	zenith=(Z0->co[I][J]-Z0->co[i][j])
                    	/sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                    if(zenith<=tan(alfa)){
                       	shadow->co[i][j]=0;
                       	y=0;
                    }else{
                      	shadow->co[i][j]=1;
                    }
               	}
         	}
     	}
  	}
}
}

/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte2(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,ii,y,I,J;
double zenith;

/*======================*/

 if(beta!=0){
   for(i=quadrata;i>0;i--){
       I=0;
       J=0;
       y=0;
       for(ii=i;ii<=quadrata;ii++){
           for(j=Z0->nch-floor(1/tan(beta)*(ii-i));
               j>=Z0->nch-floor(1/tan(beta)*(ii-i+1)) && j>0;j--){
               if(ii>(Z0->nrh+2*Z0->nch) && Z0->co[ii-(Z0->nrh+2*Z0->nch)][j]!=novalue){
                  /*shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=i}}}}}*/
                  if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==1 && I==0){
                     I=ii-(Z0->nrh+2*Z0->nch);
                     J=j;
                     y=1;
                  }else if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[ii-(Z0->nrh+2*Z0->nch)][j])
                                   /sqrt(pow((double)(I-(ii-(Z0->nrh+2*Z0->nch)))*(double)delta,(double)2)+
                                           pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=0;
                              I=ii-(Z0->nrh+2*Z0->nch);
                              J=j;
                           }else{
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=1;
                           }
                  }else if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[ii-(Z0->nrh+2*Z0->nch)][j])
                                   /sqrt(pow((double)(I-(ii-(Z0->nrh+2*Z0->nch)))*(double)delta,(double)2)+
                                           pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=0;
                              y=0;
                           }else{
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=1;
                           }
                  }
              }
        }
     }
  }
 }else{
   for(i=1;i<=Z0->nrh;i++){
       I=0;
       J=0;
       y=0;
       for(j=Z0->nch;j>0;j--){
           if(Z0->co[i][j]!=novalue){
              if(curv->co[i][j]==1 && I==0){
                     I=i;
                     J=j;
                     y=1;
              }else if(curv->co[i][j]==1 && I!=0){
                      zenith=(Z0->co[I][J]-Z0->co[i][j])
                             /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                      if(zenith<=tan(alfa)){
                         shadow->co[i][j]=0;
                         I=i;
                         J=j;
                      }else{
                         shadow->co[i][j]=1;
                      }
              }else if(curv->co[i][j]==0 && y==1){
                      zenith=(Z0->co[I][J]-Z0->co[i][j])
                              /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                      if(zenith<=tan(alfa)){
                         shadow->co[i][j]=0;
                         y=0;
                      }else{
                         shadow->co[i][j]=1;
                      }
                  }
              }
        }
     }
  }

}

/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte3(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,ii,y,I,J;
double zenith;

/*======================*/

   for(i=1;i<=quadrata;i++){
       I=0;
       J=0;
       y=0;
       for(ii=i;ii>0;ii--){
           for(j=Z0->nch-floor(1/tan(beta)*(i-ii));
               j>=Z0->nch-floor(1/tan(beta)*(i-ii+1)) && j>0;j--){
               if(ii<=Z0->nrh && Z0->co[ii][j]!=novalue){
                  /*shadow->co[ii][j]=i;}}}}*/
                  if(curv->co[ii][j]==1 && I==0){
                     I=ii;
                     J=j;
                     y=1;
                  }else if(curv->co[ii][j]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[ii][j])
                                   /sqrt(pow((double)(I-ii)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii][j]=0;
                              I=ii;
                              J=j;
                           }else{
                              shadow->co[ii][j]=1;
                           }
                  }else if(curv->co[ii][j]==0 && I!=0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[ii][j])
                                   /sqrt(pow((double)(I-ii)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii][j]=0;
                              y=0;
                           }else{
                              shadow->co[ii][j]=1;
                           }
                  }
              }
        }
     }
  }

}


/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte4(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,jj,y,I,J;
double zenith;

/*======================*/

 if(beta!=0){
   for(j=1;j<=quadrata;j++){
       I=0;
       J=0;
       y=0;
       for(jj=j;jj>0;jj--){
           for(i=Z0->nrh-floor(1/tan(beta)*(j-jj));
               i>=Z0->nrh-floor(1/tan(beta)*(j-jj+1)) && i>0;i--){
               if(jj<=Z0->nch && Z0->co[i][jj]!=novalue){
                  /*shadow->co[i][jj]=j;}}}}}*/
                  if(curv->co[i][jj]==1 && I==0){
                     I=i;
                     J=jj;
                     y=1;
                  }else if(curv->co[i][jj]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj])
                                   /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-jj)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj]=0;
                              I=i;
                              J=jj;
                           }else{
                              shadow->co[i][jj]=1;
                           }
                  }else if(curv->co[i][jj]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj])
                                   /sqrt((double)pow((I-i)*(double)delta,(double)2)+pow((double)(J-jj)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj]=0;
                              y=0;
                           }else{
                              shadow->co[i][jj]=1;
                           }
                  }
              }
        }
     }
  }
 }else{
   for(j=1;j<=Z0->nch;j++){
       I=0;
       J=0;
       y=0;
       for(i=Z0->nrh;i>0;i--){
           if(Z0->co[i][j]!=novalue){
              if(curv->co[i][j]==1 && I==0){
                 I=i;
                 J=j;
                 y=1;
              }else if(curv->co[i][j]==1 && I!=0){
                 zenith=(Z0->co[I][J]-Z0->co[i][j])
                         /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                       if(zenith<=tan(alfa)){
                          shadow->co[i][j]=0;
                          I=i;
                          J=j;
                       }else{
                          shadow->co[i][j]=1;
                        }
              }else if(curv->co[i][j]==0 && y==1){
                       zenith=(Z0->co[I][J]-Z0->co[i][j])
                               /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                       if(zenith<=tan(alfa)){
                          shadow->co[i][j]=0;
                          y=0;
                       }else{
                          shadow->co[i][j]=1;
                       }
                  }
              }
        }
     }
  }

}


/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte5(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,jj,y,I,J;
double zenith;

/*======================*/


   for(j=quadrata;j>0;j--){
       I=0;
       J=0;
       y=0;
       for(jj=j;jj<=quadrata;jj++){
           for(i=Z0->nrh-floor(1/tan(beta)*(jj-j));
               i>=Z0->nrh-floor(1/tan(beta)*(jj-j+1)) && i>0;i--){
               if(jj>quadrata-Z0->nch && Z0->co[i][jj-(quadrata-Z0->nch)]!=novalue){
                  /*shadow->co[i][jj-(quadrata-Z0->nch)]=j;}}}}*/
                  if(curv->co[i][jj-(quadrata-Z0->nch)]==1 && I==0){
                     I=i;
                     J=jj-(quadrata-Z0->nch);
                     y=1;
                  }else if(curv->co[i][jj-(quadrata-Z0->nch)]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj-(quadrata-Z0->nch)])
                                   /sqrt(pow((double)(I-i)*(double)delta,(double)2)+
                                           pow((double)(J-(jj-(quadrata-Z0->nch)))*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj-(quadrata-Z0->nch)]=0;
                              I=i;
                              J=jj-(quadrata-Z0->nch);
                           }else{
                              shadow->co[i][jj-(quadrata-Z0->nch)]=1;
                           }
                  }else if(curv->co[i][jj-(quadrata-Z0->nch)]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj-(quadrata-Z0->nch)])
                                   /sqrt(pow((double)(I-i)*(double)delta,(double)2)+
                                           pow((double)(J-(jj-(quadrata-Z0->nch)))*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj-(quadrata-Z0->nch)]=0;
                              y=0;
                           }else{
                              shadow->co[i][jj-(quadrata-Z0->nch)]=1;
                           }
                  }
              }
        }
     }
  }

}


/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte6(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,ii,y,I,J;
double zenith;

/*======================*/

 if(beta!=0){
   for(i=1;i<=quadrata;i++){
       I=0;
       J=0;
       y=0;
       for(ii=i;ii>0;ii--){
           for(j=floor(1/tan(beta)*(i-ii))+1;
               j<=floor(1/tan(beta)*(i-ii+1)) && j<=Z0->nch;j++){
               if(ii<=Z0->nrh && Z0->co[ii][j]!=novalue){
                  /*shadow->co[ii][j]=i;}}}}}*/
                  if(curv->co[ii][j]==1 && I==0){
                     I=ii;
                     J=j;
                     y=1;
                  }else if(curv->co[ii][j]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[ii][j])
                                   /sqrt(pow((double)(I-ii)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii][j]=0;
                              I=ii;
                              J=j;
                           }else{
                              shadow->co[ii][j]=1;
                           }
                  }else if(curv->co[ii][j]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[ii][j])
                                   /sqrt((double)pow((I-ii)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii][j]=0;
                              y=0;
                           }else{
                              shadow->co[ii][j]=1;
                           }
                  }
              }
        }
     }
  }
 }else{
   for(i=1;i<=Z0->nrh;i++){
       I=0;
       J=0;
       y=0;
       for(j=1;j<=Z0->nch;j++){
           if(Z0->co[i][j]!=novalue){
              if(curv->co[i][j]==1 && I==0){
                 I=i;
                 J=j;
                 y=1;
              }else if(curv->co[i][j]==1 && I!=0){
                 zenith=(Z0->co[I][J]-Z0->co[i][j])
                        /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                 if(zenith<=tan(alfa)){
                    shadow->co[i][j]=0;
                    I=i;
                    J=j;
                 }else{
                    shadow->co[i][j]=1;
                    }
              }else if(curv->co[i][j]==0 && y==1){
                    zenith=(Z0->co[I][J]-Z0->co[i][j])
                           /sqrt(pow((double)(I-i)*(double)delta,(double)2)+pow((double)(J-j)*(double)delta,(double)2));
                    if(zenith<=tan(alfa)){
                       shadow->co[i][j]=0;
                       y=0;
                    }else{
                       shadow->co[i][j]=1;
                    }
                  }
              }
        }
     }
  }

}


/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte7(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,ii,y,I,J;
double zenith;

/*======================*/


   for(i=quadrata;i>0;i--){
       I=0;
       J=0;
       y=0;
       for(ii=i;ii<quadrata;ii++){
           for(j=floor(1/tan(beta)*(ii-i))+1;
               j<=floor(1/tan(beta)*(ii-i+1)) && j<=Z0->nch;j++){
               if(ii>(Z0->nrh+2*Z0->nch) && Z0->co[ii-(Z0->nrh+2*Z0->nch)][j]!=novalue){
                  /*shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=i;}}}}*/
                  if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==1 && I==0){
                     I=ii-(Z0->nrh+2*Z0->nch);
                     J=j;
                     y=1;
                  }else if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==1 && I!=0){
                           zenith=(Z0->co[I][J]-Z0->co[ii-(Z0->nrh+2*Z0->nch)][j])
                                   /sqrt(pow((double)(I-(ii-(Z0->nrh+2*Z0->nch)))*(double)delta,(double)2)+
                                           pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=0;
                              I=ii-(Z0->nrh+2*Z0->nch);
                              J=j;
                           }else{
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=1;
                           }
                  }else if(curv->co[ii-(Z0->nrh+2*Z0->nch)][j]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[ii-(Z0->nrh+2*Z0->nch)][j])
                                   /sqrt(pow((double)(I-(ii-(Z0->nrh+2*Z0->nch)))*(double)delta,(double)2)+
                                           pow((double)(J-j)*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=0;
                              y=0;
                           }else{
                              shadow->co[ii-(Z0->nrh+2*Z0->nch)][j]=1;
                           }
                  }
              }
        }
     }
  }

}


/*----------------------------------------------------------------------------------------------------------*/
void Orizzonte8(double delta,long quadrata,double beta,double alfa,DOUBLEMATRIX *Z0,
                SHORTMATRIX *curv,SHORTMATRIX *shadow, double novalue)
{
/*=====================*/

long i,j,jj,y,I,J;
double zenith;

/*======================*/


   for(j=quadrata;j>0;j--){
       I=0;
       J=0;
       y=0;
       for(jj=j;jj<=quadrata;jj++){
           for(i=floor(1/tan(beta)*(jj-j))+1;
               i<=floor(1/tan(beta)*(jj-j+1)) && i<=Z0->nrh;i++){
               if(jj>quadrata-Z0->nch && Z0->co[i][jj-(quadrata-Z0->nch)]!=novalue){
                  /*shadow->co[i][jj-(quadrata-Z0->nch)]=j;}}}}*/
                  if(curv->co[i][jj]==1 && I==0){
                     I=i;
                     J=jj-(quadrata-Z0->nch);
                     y=1;
                  }else if(curv->co[i][jj-(quadrata-Z0->nch)]==1 && I!=1){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj-(quadrata-Z0->nch)])
                                   /sqrt(pow((double)(I-i)*(double)delta,(double)2)+
                                           pow((double)(J-(jj-(quadrata-Z0->nch)))*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj-(quadrata-Z0->nch)]=0;
                              I=i;
                              J=jj-(quadrata-Z0->nch);
                           }else{
                              shadow->co[i][jj-(quadrata-Z0->nch)]=1;
                           }
                  }else if(curv->co[i][jj-(quadrata-Z0->nch)]==0 && y==1){
                           zenith=(Z0->co[I][J]-Z0->co[i][jj-(quadrata-Z0->nch)])
                                   /sqrt(pow((double)(I-i)*(double)delta,(double)2)+
                                           pow((double)(J-(jj-(quadrata-Z0->nch)))*(double)delta,(double)2));
                           if(zenith<=tan(alfa)){
                              shadow->co[i][jj-(quadrata-Z0->nch)]=0;
                              y=0;
                           }else{
                              shadow->co[i][jj-(quadrata-Z0->nch)]=1;
                           }
                  }
              }
        }
     }
  }

}



