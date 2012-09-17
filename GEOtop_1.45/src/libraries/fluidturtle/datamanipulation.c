//#include "turtle.h"
#include "t_datamanipulation.h"
//#include "t_alloc.h"
//#include "t_statistics.h"
//#include "write_dem.h"

/*-------------------------------------------------------------------------------*/



DOUBLEVECTOR *vectorize_doublematrix(DOUBLEMATRIX *input)

{

  DOUBLEVECTOR *U;

  U=(DOUBLEVECTOR *)malloc(sizeof(DOUBLEVECTOR));
  if(!U){
    t_error("This vector cannot be allocated");
  }

  U->nl=1;
  U->nh=(input->nrh)*(input->nch);
  U->isdynamic=1;
  if(input->co==NULL || input->isdynamic!=1 ){

    t_error("This matrix is not properly allocated");

  } else {

    U->co=input->co[1];

  }

  /* print_doublevector_elements(U,10); */


  free(input);

  return U;


}


/*-------------------------------------------------------------------------------*/



LONGVECTOR *vectorize_longmatrix(LONGMATRIX *input)

{

  LONGVECTOR *U;

  U=(LONGVECTOR *)malloc(sizeof(LONGVECTOR));
  if(!U){
    t_error("This vector cannot be allocated");
  }

  U->nl=1;
  U->nh=(input->nrh)*(input->nch);
  U->isdynamic=1;
  if(input->co==NULL || input->isdynamic!=1 ){

    t_error("This matrix is not properly allocated");

  } else {

    U->co=input->co[1];

  }

  /* print_doublevector_elements(U,10); */


  free(input);

  return U;


}
/*-------------------------------------------------------------------------------*/



FLOATVECTOR *vectorize_floatmatrix(FLOATMATRIX *input)

{

  FLOATVECTOR *U;

  U=(FLOATVECTOR *)malloc(sizeof(FLOATVECTOR));
  if(!U){
    t_error("This vector cannot be allocated");
  }

  U->nl=1;
  U->nh=(input->nrh)*(input->nch);
  U->isdynamic=1;
  if(input->co==NULL || input->isdynamic!=1 ){

    t_error("This matrix is not properly allocated");

  } else {

    U->co=input->co[1];

  }

  /* print_floatvector_elements(U,10); */


  free(input);

  return U;


}


/*-------------------------------------------------------------------------------*/
SHORTVECTOR *vectorize_shortmatrix(SHORTMATRIX *input)

{

  SHORTVECTOR *U;

  U=(SHORTVECTOR *)malloc(sizeof(SHORTVECTOR));
  if(!U){
    t_error("This vector cannot be allocated");
  }

  U->nl=1;
  U->nh=(input->nrh)*(input->nch);
  U->isdynamic=1;
  if(input->co==NULL || input->isdynamic!=1 ){

    t_error("This matrix is not properly allocated");

  } else {

    U->co=input->co[1];

  }

  /* print_floatvector_elements(U,10); */


  free(input);

  return U;


}


/*-------------------------------------------------------------------------------*/

void sortreal(DOUBLEVECTOR *ra)

{

  long n;
  unsigned int i,ir,j,l;
  double rra;

  n=ra->nh;

  if (n <2) return;
  l=(n>>1)+1;
  ir=n;

  for(;;){
    if (l >1) {
      rra=ra->co[--l];
    } else {
      rra=ra->co[ir];
      ra->co[ir]=ra->co[1];
      if (--ir ==1){
	ra->co[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j<=ir){
      if(j<ir && ra->co[j] < ra->co[j+1]) j++;
      if(rra < ra->co[j]){
	ra->co[i]=ra->co[j];
	i=j;
	j <<=1;
      } else j=ir+1;
    }
    ra->co[i]=rra;
  }

}

/*-------------------------------------------------------------------------------*/
void  sort2realvectors(DOUBLEVECTOR *ra,DOUBLEVECTOR *rb)

{

  long n;
  unsigned int i,ir,j,l;
  double rra,rrb;

  n=ra->nh;

  if (n <2) return;
  l=(n>>1)+1;
  ir=n;

  for(;;){
    if (l >1) {
      rra=ra->co[--l];
      rrb=rb->co[l];
    } else {
      rra=ra->co[ir];
      ra->co[ir]=ra->co[1];
      rrb=rb->co[ir];
      rb->co[ir]=rb->co[1];
	 
      if (--ir ==1){
	ra->co[1]=rra;
	rb->co[1]=rrb;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j<=ir){
      if(j<ir && ra->co[j] < ra->co[j+1]) j++;
      if(rra < ra->co[j]){
	ra->co[i]=ra->co[j];
	rb->co[i]=rb->co[j];
	i=j;
	j <<=1;
      } else j=ir+1;
    }
    ra->co[i]=rra;
    rb->co[i]=rrb;
    
  }


}

/*--------------------------------------------------------------------------------*/

void  sort2floatvectors(FLOATVECTOR *ra,FLOATVECTOR *rb)

{

  long n;
  unsigned int i,ir,j,l;
  float  rra;
  float rrb;

  n=ra->nh;

  if (n <2) return;
  l=(n>>1)+1;
  ir=n;
  for(;;){
    if (l >1) {
      rra=ra->co[--l];
      rrb=rb->co[l];
    } else {
      rra=ra->co[ir];
      ra->co[ir]=ra->co[1];
      rrb=rb->co[ir];
      rb->co[ir]=rb->co[1];
	 
      if (--ir ==1){
	ra->co[1]=rra;
	rb->co[1]=rrb;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j<=ir){
      if(j<ir && ra->co[j] < ra->co[j+1]) j++;
      if(rra < ra->co[j]){
	ra->co[i]=ra->co[j];
	rb->co[i]=rb->co[j];
	i=j;
	j <<=1;
      } else j=ir+1;
    }
    ra->co[i]=rra;
    rb->co[i]=rrb;
    
  }

}

/*-------------------------------------------------------------------------------*/
void  sort2vectors(LONGVECTOR *ra,FLOATVECTOR *rb)

{

  long n;
  unsigned int i,ir,j,l;
  long rra;
  float rrb;

  n=ra->nh;

  if (n <2) return;
  l=(n>>1)+1;
  ir=n;

  for(;;){
    if (l >1) {
      rra=ra->co[--l];
      rrb=rb->co[l];
    } else {
      rra=ra->co[ir];
      ra->co[ir]=ra->co[1];
      rrb=rb->co[ir];
      rb->co[ir]=rb->co[1];
	 
      if (--ir ==1){
	ra->co[1]=rra;
	rb->co[1]=rrb;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j<=ir){
      if(j<ir && ra->co[j] < ra->co[j+1]) j++;
      if(rra < ra->co[j]){
	ra->co[i]=ra->co[j];
	rb->co[i]=rb->co[j];
	i=j;
	j <<=1;
      } else j=ir+1;
    }
    ra->co[i]=rra;
    rb->co[i]=rrb;
    
  }


}



/*-------------------------------------------------------------------------------*/
/*
Inputs: 1) a sorted vector containing the data to be binned; 2) the number of bins
(instead of fixing the size of each bin is usually convenient to select a fixed number 
of bins ). If the number of bins is set to 0, is simply counted the number of elements with the
same abscissa; 3) the minimum number of elements required in each bin (this implies
that if a bin has not enough elements the numbers of elements of two adjacent bins are summed).
In the case of the exponential binning, the j-th bin extends from base^(j*delta) 
to base^((j+1)*delta) where base is the fourth input field and delta the extension of 
each bin in the logarithm  axis.


Return: a matrix of double: the first column contains the number of elements
in each mean, the second coulumns the mean abscissa of the data in the bin, the third
the highest limit of each bin interval */

DOUBLEMATRIX *simplehystogram(DOUBLEVECTOR *U,long N,long mn)

{

  long j,count,count1;
  double delta,min,max,tmp;
  DOUBLEMATRIX *indx=NULL;
  XYZ *bins,*head,*pre;
  /* char ch; */

  if(N <0 ) t_error("A negative number of bins is not allowed");

    
  head=new_xyz();
  head->y=U->co[1];
  pre=head;
  bins=head;
  count1=1;
  count=2;
  while(count <=U->nh){
       
    while(U->co[count]==U->co[count-1] && count <=U->nh){
      count++;
    }
		
    bins->next=new_xyz();
    bins=bins->next;
    bins->x=count-count1;
    count1=count;
    head->x++;
    bins->y=U->co[count-1];
    count++;
	
  }
	

  bins=head->next;

  while(bins!=NULL){
    count=bins->x;
    while(bins->x < mn && bins->next!=NULL){
      bins->y=(bins->y*bins->x+bins->next->y*bins->next->x)/(bins->x+bins->next->x);
      bins->x+=bins->next->x;
      bins->z=0;
      delete_xyz(head,bins);
      (head->x)--;
    }
	
    if(bins->x< mn) {
	
      if(bins!=head){
	pre->x+=bins->x;
	delete_xyz(head,pre);
      }else {
	printf("\nWarning::Obtaining just one bin\n");
      }
	
    }
	
    pre=bins;
    bins=bins->next;
  }

  /* print_xyz_elements(head,10); */

  if(N<1){

    indx=new_doublematrix(head->x,3);
    xyz_into_doublematrix(head,indx);
    delete_xyz_list(head);

  }else{

    if(head->x <N) printf("\nWarning::Some bin are remaining empty\n");
    /*	cnt=count_xyz_elements(head); */
    indx=new_doublematrix(N,3);
    delta=( U->co[U->nh]-head->next->y)/(N-1);
    max=U->co[U->nh]+0.5*delta;
    min=head->next->y-0.5*delta;
    head->y=head->z=min;
    bins=head->next;    
    for(j=1;j<=N;j++){
      indx->co[j][1]=0;
      indx->co[j][2]=0;
      indx->co[j][3]=min+j*delta;
      tmp=0;	
      if(bins!=NULL){	
	while(bins->y<=min+j*delta){
	  tmp+=bins->x;
	  indx->co[j][2]=(indx->co[j][1]*indx->co[j][2]+bins->y*bins->x)/tmp;
	  indx->co[j][1]=tmp;
	  bins=bins->next;
	  if(bins==NULL) break;
	}		
      }

    }

    delete_xyz_list(head);	

  }

  /* print_doublematrix_elements(indx,10); */

  return indx;

}

/*----------------------------------------------------------------------*/

DOUBLEMATRIX *exponentialhystogram(DOUBLEVECTOR *U,long N,long mn,long base)

{

  long j,count,count1;
  double logdelta,min,max,tmp,logbase;
  DOUBLEMATRIX *indx=NULL;
  XYZ *bins,*head,*pre;


  if(N <0 ) t_error("A negative number of bins is not allowed");

  logbase=log(base);
  head=new_xyz();
  head->y=U->co[1];
  pre=head;
  bins=head;
  count1=1;
  count=2;
  while(count <=U->nh){
       
    while(U->co[count]==U->co[count-1] && count <=U->nh){
      count++;
    }
		
    bins->next=new_xyz();
    bins=bins->next;
    bins->x=count-count1;
    count1=count;
    head->x++;
    bins->y=U->co[count-1];
    count++;
	
  }
	

  bins=head->next;

  while(bins!=NULL){
    count=bins->x;
    while(bins->x < mn && bins->next!=NULL){
      bins->y=(bins->y*bins->x+bins->next->y*bins->next->x)/(bins->x+bins->next->x);
      bins->x+=bins->next->x;
      bins->z=0;
      delete_xyz(head,bins);
      (head->x)--;
    }
	
    if(bins->x< mn) {
	
      if(bins!=head){
	pre->x+=bins->x;
	delete_xyz(head,pre);
      }else {
	printf("\nWarning::Obtaining just one bin\n");
      }
	
    }
	
    pre=bins;
    bins=bins->next;
  }

  /* print_xyz_elements(head,10); */


  if(N<1){

    indx=new_doublematrix(head->x,3);
    xyz_into_doublematrix(head,indx);
    /*	print_doublematrix_elements(indx,10); */
    delete_xyz_list(head);

  }else{

    if(head->x <N) printf("\nWarning::Some bin is remaining empty\n");
    /*	cnt=count_xyz_elements(head); */
    indx=new_doublematrix(N,3);
    logdelta=(log( U->co[U->nh])-log(head->next->y))/(logbase*(N-1));
    max=log(U->co[U->nh])/log(base)+0.5*logdelta;
    min=log(head->next->y)/logbase-0.5*logdelta;
    head->y=head->z=min;
    bins=head->next;    

    for(j=1;j<=N;j++){
      indx->co[j][1]=0;
      indx->co[j][2]=0;
      indx->co[j][3]=pow(base,min+j*logdelta);
      tmp=0;	
      if(bins!=NULL){	
	while(log(bins->y)/logbase<=min+j*logdelta){
	  tmp+=bins->x;
	  indx->co[j][2]=(indx->co[j][1]*indx->co[j][2]+bins->y*bins->x)/tmp;
	  indx->co[j][1]=tmp;
	  bins=bins->next;
	  if(bins==NULL) break;
	}		
      }
    }

    delete_xyz_list(head);	

  }

  /* print_doublematrix_elements(indx,10); */

  return indx;

}


/*----------------------------------------------------------------------*/
void realpair_into_doublematrix(REALPAIR * head,DOUBLEMATRIX *indx )
{

  long i;
  REALPAIR *tmp;

  if(head==NULL || indx==NULL){

    t_error("Inputs are not properly allocated");
	
  } else if(head->x!=indx->nrh) {

    t_error("Inputs do not have the same dimensions");
  }

  tmp=head->next;

  for(i=1;i<=head->x;i++){
    if(tmp==NULL) t_error("Something was wrong here");
    indx->co[i][1]=tmp->x;
    indx->co[i][2]=tmp->y;
    tmp=tmp->next;
  }

}

/*----------------------------------------------------------------------*/


void xyz_into_doublematrix(XYZ * head,DOUBLEMATRIX *indx )
{

  long i;
  XYZ *tmp;

  if(head==NULL || indx==NULL){

    t_error("Inputs are not properly allocated");
	
  } else if(head->x!=indx->nrh) {

    t_error("Inputs do not have the same dimensions");
  }

  tmp=head->next;

  for(i=1;i<=head->x;i++){
    if(tmp==NULL) t_error("Something was wrong here");
    indx->co[i][1]=tmp->x;
    indx->co[i][2]=tmp->y;
    indx->co[i][3]=tmp->z;
    tmp=tmp->next;
  }

}

/*---------------------------------------------------------------------------*/
void initialize_longvector(LONGVECTOR *L,long sign)

{

  long i;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nh;i++){
	L->co[i]=sign;
      }
    }else{
      t_error("This longvector was no properly allocated");
    }
  }else {
    t_error("A null vector was addressed");
  }
}

/*---------------------------------------------------------------------------*/
void initialize_shortvector(SHORTVECTOR *L,short sign)

{

  long i;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nh;i++){
	L->co[i]=sign;
      }
    }else{
      t_error("This shortvector was no properly allocated");
    }
  }else {
    t_error("A null vector was addressed");
  }
}

/*---------------------------------------------------------------------------*/
void initialize_doublevector(DOUBLEVECTOR *L,double sign)

{

  long i;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=L->nl;i<=L->nh;i++){
	L->co[i]=sign;
      }
    }else{
      t_error("This doublevector was no properly allocated");
    }
  }else {
    t_error("A null vector was addressed");
  }
}

/*---------------------------------------------------------------------------*/
void initialize_floatvector(FLOATVECTOR *L,float sign)

{

  long i;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nh;i++){
	L->co[i]=sign;
      }
    }else{
      t_error("This floatvector was no properly allocated");
    }
  }else {
    t_error("A null vector was addressed");
  }
}

/*---------------------------------------------------------------------------*/
void initialize_shortmatrix(SHORTMATRIX *L,short sign)

{

  long i,j;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nrh;i++){
	for(j=1;j<=L->nch;j++){
	  L->co[i][j]=sign;
	}
      }
    }else{
      t_error("This matrix was no properly allocated");
    }
  }else {
    t_error("A null matrix was addressed");
  }
}


/*---------------------------------------------------------------------------*/
void initialize_longmatrix(LONGMATRIX *L,long sign)

{

  long i,j;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nrh;i++){
	for(j=1;j<=L->nch;j++){
	  L->co[i][j]=sign;
	}
      }
    }else{
      t_error("This matrix was no properly allocated");
    }
  }else {
    t_error("A null matrix was addressed");
  }
}


/*---------------------------------------------------------------------------*/
void initialize_floatmatrix(FLOATMATRIX *L,float sign)

{

  long i,j;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nrh;i++){
	for(j=1;j<=L->nch;j++){
	  L->co[i][j]=sign;
	}
      }
    }else{
      t_error("This matrix was no properly allocated");
    }
  }else {
    t_error("A null matrix was addressed");
  }
}

/*---------------------------------------------------------------------------*/
void initialize_doublematrix(DOUBLEMATRIX *L,double sign)

{

  long i,j;

  if(L!=NULL){
    if(L->isdynamic==1){
      for(i=1;i<=L->nrh;i++){
	for(j=1;j<=L->nch;j++){
	  L->co[i][j]=sign;
	}
      }
    }else{
      t_error("This matrix was no properly allocated");
    }
  }else {
    t_error("A null matrix was addressed");
  }
}


/*--------------------------------------------------------------------------*/

void copy_shortmatrix(SHORTMATRIX *origin,SHORTMATRIX *destination)

{

  long i,j;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A matrix was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){

    t_error("A matrix was not allocated properly");

  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ){
	
    t_error("The matrixes do not have the same dimensions");

  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
	
      destination->co[i][j]=origin->co[i][j];

    }
  }

}

/*--------------------------------------------------------------------------*/

void copy_intmatrix(INTMATRIX *origin,INTMATRIX *destination)

{

  long i,j;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A matrix was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){

    t_error("A matrix was not allocated properly");

  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ){
	
    t_error("The matrixes do not have the same dimensions");

  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
	
      destination->co[i][j]=origin->co[i][j];

    }
  }

}





/*--------------------------------------------------------------------------*/

void copy_longmatrix(LONGMATRIX *origin,LONGMATRIX *destination)

{

  long i,j;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A matrix was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){

    t_error("A matrix was not allocated properly");

  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ){
	
    t_error("The matrixes do not have the same dimensions");

  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
	
      destination->co[i][j]=origin->co[i][j];

    }
  }

}



/*--------------------------------------------------------------------------*/

void copy_floatmatrix(FLOATMATRIX *origin,FLOATMATRIX *destination)

{

  long i,j;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A matrix was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){

    t_error("A matrix was not allocated properly");

  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ){
	
    t_error("The matrixes do not have the same dimensions");

  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
	
      destination->co[i][j]=origin->co[i][j];

    }
  }

}


/*--------------------------------------------------------------------------*/

void copy_doublematrix(DOUBLEMATRIX *origin,DOUBLEMATRIX *destination)

{

  long i,j;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A matrix was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nrh <1 || destination->nrh <1 ||  origin->nch <1 || destination->nch <1 ){

    t_error("A matrix was not allocated properly");

  }else if( origin->nrh != destination->nrh ||  origin->nch != destination->nch ){
	
    t_error("The matrixes do not have the same dimensions");

  }

  for(i=1;i<=origin->nrh;i++){
    for(j=1;j<=origin->nch;j++){
	
      destination->co[i][j]=origin->co[i][j];

    }
  }

}




/*--------------------------------------------------------------------------*/

void copy_shortvector(SHORTVECTOR *origin,SHORTVECTOR *destination)

{

  long i;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A vector was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nh <1 || destination->nh <1  ){

    t_error("A vector was not allocated properly");

  }else if( origin->nh != destination->nh  ){
	
    t_error("The vector do not have the same dimensions");

  }

  for(i=1;i<=origin->nh;i++){
	
	
    destination->co[i]=origin->co[i];

	
  }

}

/*--------------------------------------------------------------------------*/

void copy_intvector(INTVECTOR *origin,INTVECTOR *destination)

{

  long i;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A vector was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nh <1 || destination->nh <1  ){

    t_error("A vector was not allocated properly");

  }else if( origin->nh != destination->nh  ){
	
    t_error("The vector do not have the same dimensions");

  }

  for(i=1;i<=origin->nh;i++){
	
	
    destination->co[i]=origin->co[i];

	
  }

}


/*--------------------------------------------------------------------------*/

void copy_longvector(LONGVECTOR *origin,LONGVECTOR *destination)

{

  long i;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A vector was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nh <1 || destination->nh <1 ){

    t_error("A vector was not allocated properly");

  }else if( origin->nh != destination->nh  ){
	
    t_error("The vector do not have the same dimensions");

  }

  for(i=1;i<=origin->nh;i++){
	
	
    destination->co[i]=origin->co[i];

	
  }

}


/*--------------------------------------------------------------------------*/

void copy_floatvector(FLOATVECTOR *origin,FLOATVECTOR *destination)

{

  long i;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A vector was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nh <1 || destination->nh <1 ){

    t_error("A vector was not allocated properly");

  }else if( origin->nh != destination->nh  ){
	
    t_error("The vector do not have the same dimensions");

  }

  for(i=1;i<=origin->nh;i++){
	
	
    destination->co[i]=origin->co[i];

	
  }

}

/*--------------------------------------------------------------------------*/

void copy_doublevector(DOUBLEVECTOR *origin,DOUBLEVECTOR *destination)

{

  long i;

  if(origin==NULL || destination==NULL || origin->co==NULL || destination->co==NULL){

    t_error("A vector was not allocated");

  } else if(origin->isdynamic!=1 || destination->isdynamic!=1 || origin->nh <1 || destination->nh <1  ){

    t_error("A vector was not allocated properly");

  }else if( origin->nh != destination->nh  ){
	
    t_error("The vector do not have the same dimensions");

  }

  for(i=1;i<=origin->nh;i++){
	
	
    destination->co[i]=origin->co[i];

	
  }

}

/*--------------------------------------------------------------------------*/

void add_doublevector(DOUBLEVECTOR *small, DOUBLEVECTOR *big)

{
	
	long i;
	
	if(small==NULL || big==NULL || small->co==NULL || big->co==NULL){
		
		t_error("A vector was not allocated");
		
	} else if(small->isdynamic!=1 || big->isdynamic!=1 || small->nh <1 || big->nh <1  ){
		
		t_error("A vector was not allocated properly");
		
	}else if( small->nh != big->nh  ){
		
		t_error("The vector do not have the same dimensions");
		
	}
	
	for(i=1;i<=small->nh;i++){
		
		
		big->co[i] += small->co[i];
		
		
	}
	
}

/*--------------------------------------------------------------------------*/
void floatmatrix_element_multiplication(FLOATMATRIX* dist,FLOATMATRIX *ca)

{

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
      dist->co[i][j]*=ca->co[i][j]; 	
    }
  }


}

/*--------------------------------------------------------------------------*/
void doublematrix_element_multiplication(DOUBLEMATRIX* dist,DOUBLEMATRIX *ca)

{

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
      dist->co[i][j]*=ca->co[i][j]; 	
    }
  }


}

/*--------------------------------------------------------------------------*/
void longmatrix_element_multiplication(LONGMATRIX* dist,LONGMATRIX *ca)

{

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
      dist->co[i][j]*=ca->co[i][j]; 	
    }
  }


}


/*--------------------------------------------------------------------------*/
void shortmatrix_element_multiplication(SHORTMATRIX* dist,SHORTMATRIX *ca)

{

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
      dist->co[i][j]*=ca->co[i][j]; 	
    }
  }


}

/*--------------------------------------------------------------------------*/
void shortvector_element_multiplication(SHORTVECTOR* dist,SHORTVECTOR *ca)

{

  long i;

  if(dist==NULL || ca==NULL || dist->co==NULL || ca->co==NULL){
    t_error("Vectors were not allocated");
  } else if(dist->isdynamic!=1 || dist->isdynamic!=1 || dist->nh <2 ){
    t_error("Vectors were not properly allocated");
  } else if(dist->nh !=ca->nh ) {
    t_error("These vectors do not have the same dimensions");
  }

  for(i=1;i<=dist->nh;i++){
    dist->co[i]*=ca->co[i]; 	
  }


}

/*--------------------------------------------------------------------------*/
void longvector_element_multiplication(LONGVECTOR* dist,LONGVECTOR *ca)

{

  long i;

  if(dist==NULL || ca==NULL || dist->co==NULL || ca->co==NULL){
    t_error("Vectors were not allocated");
  } else if(dist->isdynamic!=1 || dist->isdynamic!=1 || dist->nh <2 ){
    t_error("Vectors were not properly allocated");
  } else if(dist->nh !=ca->nh ) {
    t_error("These vectors do not have the same dimensions");
  }

  for(i=1;i<=dist->nh;i++){
    dist->co[i]*=ca->co[i]; 	
  }


}

/*--------------------------------------------------------------------------*/
void floatvector_element_multiplication(FLOATVECTOR* dist,FLOATVECTOR *ca)

{

  long i;

  if(dist==NULL || ca==NULL || dist->co==NULL || ca->co==NULL){
    t_error("Vectors were not allocated");
  } else if(dist->isdynamic!=1 || dist->isdynamic!=1 || dist->nh <2 ){
    t_error("Vectors were not properly allocated");
  } else if(dist->nh !=ca->nh ) {
    t_error("These vectors do not have the same dimensions");
  }

  for(i=1;i<=dist->nh;i++){
    dist->co[i]*=ca->co[i]; 	
  }


}

/*--------------------------------------------------------------------------*/
void doublevector_element_multiplication(DOUBLEVECTOR* dist,DOUBLEVECTOR *ca)

{

  long i;

  if(dist==NULL || ca==NULL || dist->co==NULL || ca->co==NULL){
    t_error("Vectors were not allocated");
  } else if(dist->isdynamic!=1 || dist->isdynamic!=1 || dist->nh <2 ){
    t_error("Vectors were not properly allocated");
  } else if(dist->nh !=ca->nh ) {
    t_error("These vectors do not have the same dimensions");
  }

  for(i=1;i<=dist->nh;i++){
    dist->co[i]*=ca->co[i]; 	
  }


}
/*-----------------------------------------------------------------------*/


void clean_floatmatrix(FLOATMATRIX *iv,FLOATMATRIX *ov,FLOATVECTOR *U,FLOATVECTOR *V)


{

  long i,j;

  if(ov->nrh!=iv->nrh || ov->nch!=iv->nch){

    t_error("Matrixes do not have proper dimensions");

  }

  if(V->co[1] >0){
    for(i=1;i<=iv->nrh;i++){
      for(j=1;j<=iv->nch;j++){
	if(ov->co[i][j] >= V->co[2]){
	  iv->co[i][j]=U->co[2];
	}
       		
					
      }
    }
  }else if(V->co[1] <0){

    for(i=1;i<=iv->nrh;i++){
      for(j=1;j<=iv->nch;j++){
	if(ov->co[i][j] <= V->co[2]){
	  iv->co[i][j]=U->co[2];
	}
					
      }
    }

  }else{

    for(i=1;i<=iv->nrh;i++){
      for(j=1;j<=iv->nch;j++){
	if(ov->co[i][j] == V->co[2]){
	  iv->co[i][j]=U->co[2];
	}
					
      }
    }


  }

}

/*---------------------------------------------------------------------------*/

DOUBLEBIN *split(DOUBLEVECTOR *U,long N,FLOATVECTOR *novalue)


{

  double delta,min,max;
  long i,count,count1,minposition,maxposition,bin_vuoti,num_max;
  DOUBLEBIN *l;
  LONGPAIR *head,*bins;
  /* char ch; */

  minposition=1;
  maxposition=U->nh;

  if(N <=1){


    head=new_longpair();
    head->i=0;
    bins=head;
    count1=1;
    count=2;
    printf("This options is memory expensive\n");
    printf("ENTER THE MAXIMUM NUMBER OF ALLOWED BIN\n");
    scanf("%ld",&num_max);

    while(count <= U->nh){
      while(U->co[count]==U->co[count-1] && count <=U->nh){
	count++;
      }
      /*		printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
			scanf("%c",&ch);
      */
	
      bins->next=new_longpair();
      bins=bins->next;
      bins->i=count-count1;
      bins->j=1;
      count1=count-1;
      head->i++;
      count1=count;
      count++; 
      if(head->i>num_max)t_error("Il numero di bin eccede il numero massimo che hai consentito");

    }

  }else if(N >1){
    if(novalue->co[1]==0){
      minposition=1;
      maxposition=U->nh;
    }else if(novalue->co[1]==-1){
      minposition=2;
      /*	    printf("*****%d\n",novalue); */
      max=U->co[U->nh];
      while(U->co[minposition]<=novalue->co[2]){
	minposition++;
      }
      min=U->co[minposition];
      /*printf("***%f %f\n",min,U->co[minposition-1]);
	scanf(" %c",&ch);*/
      maxposition=U->nh;
    }else{
      maxposition=U->nh-1;
      while(U->co[maxposition]>=novalue->co[2]){
	/*			printf("%d\n",maxposition);
				scanf("%c",&ch);	
	*/			maxposition--;
      }	    
      max=U->co[maxposition];
      minposition=1;
      min=U->co[1];
    }

    printf("THE MINIMUM VALUE IS %f\n",min);
    printf("THE MAXIMUM VALUE IS %f\n",max);
    delta=(max-min)/(N-1);
    printf("DELTA IS%f\n",delta);
    /*  scanf("%c",&ch); */
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=minposition;
    count=minposition+1;
    bin_vuoti=0;

    while(count <= maxposition){
      if(U->co[count] < min+0.5*delta){
	while(U->co[count] < min+0.5*delta  && count<=maxposition){
	  count++;
	}

	bins->next=new_longpair();
	bins=bins->next;	
	bins->i=count-count1;
	/*			printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
				scanf("%c",&ch);
	*/
	bins->j=1;
	count1=count-1;
	(head->i)++;
	count1=count;
	count++; 

      }else{
	bin_vuoti++;
	printf("\nWarning:: An empty bin was encountered\n");
      }			
      min+=delta;
    }		
    if(bin_vuoti!=0){
      printf("\nThere was look for (%ld) empty bins\n",bin_vuoti);
    }
  }

  l=(DOUBLEBIN *)malloc(sizeof(DOUBLEBIN));
  if (!l) t_error("allocation failure in new_doublebin()");
  l->isdynamic=isDynamic;
  if(head->i < 2 ){
    t_error("Something wrong happened in binning");
  }else{
    l->index=new_longvector(head->i);
    (l->index)->nl=1;
    (l->index)->nh=head->i;
    bins=head->next;
    for(i=1;i<=head->i;i++){
      (l->index)->co[i]=bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }
    l->co=(double **)malloc((size_t) ((l->index->nh+NR_END)*sizeof(double *)));
    if (!l->co) t_error("allocation failure in new_doublebin()");
    l+=NR_END-1;

    l->co[1]=U->co+minposition-1;
    if(!(l->co[1])) t_error("Failure in splitting this vector");
    l->co[1]+=NR_END;
    l->co[1]-=1;
    bins=head->next;	
    for(i=2;i<=l->index->nh; i++){
      l->co[i]=l->co[i-1]+bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }

  }


  U->isdynamic=-1;
  U->nl=-1;
  U->nh=1;

  free(U);

  delete_longpair_list(head);

  return l;


}


/*---------------------------------------------------------------------------*/

DOUBLEBIN *esponentialsplit(DOUBLEVECTOR *U,long N,double base,FLOATVECTOR* novalue)


{

  double min,max,logdelta,logbase,logmin;
  long i,count,count1,minposition,maxposition,bin_vuoti,num_max;
  DOUBLEBIN *l;
  LONGPAIR *head,*bins;
  /*char ch;

  minposition=1;
  maxposition=U->nh;
  printf("%d\n",N);
  scanf("%c",&ch);*/
  logbase=log10(base);

  if(N <=1){
    printf("This options is memory expensive\n");
    printf("ENTER THE MAXIMUM NUMBER OF ALLOWED BIN\n");
    scanf("%ld",&num_max);
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=1;
    count=2;

    while(count <= U->nh){
      while(U->co[count]==U->co[count-1] && count <=U->nh){
	count++;
      }
	
      bins->next=new_longpair();
      bins=bins->next;
      bins->i=count-count1;
      /*		printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
			scanf("%c",&ch);
      */
      bins->j=1;
      count1=count-1;
      head->i++;
      count1=count;
      count++; 
      if(head->i>num_max)t_error("Il numero di bin eccede il numero massimo che hai consentito");
    }
  }else if(N >1){
    if(novalue->co[1]==0){
      minposition=1;
      if(U->co[minposition] <=0){
	while(U->co[minposition]<=0){
	  minposition++;
	  if( minposition > U->nh) t_error("No valid data for esponential binning");
	}	    
      }
      maxposition=U->nh;
    }else if(novalue->co[1]==-1){
      minposition=2;
      /*	    printf("*****%d\n",novalue); */
      max=U->co[U->nh];
      while(U->co[minposition]<=novalue->co[2]){
	minposition++;
      }
        
      min=U->co[minposition];
      maxposition=U->nh;
    }else{
      maxposition=U->nh-1;
      while(U->co[maxposition]>=novalue->co[2]){
	/*			printf("%d\n",maxposition);
				scanf("%c",&ch);	
	*/			maxposition--;
      }	    
      max=U->co[maxposition];
      minposition=1;
      min=U->co[1];
    }

    printf("THE MINIMUM VALUE IS %f\n",min);
    printf("THE MAXIMUM VALUE IS %f\n",max);
    logdelta=(log10(max)-log10(min))/(logbase*(N-1));
    printf("LOG_10  DELTA IS%f\n",logdelta);
    /*  scanf("%c",&ch); */
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=minposition;
    count=minposition+1;
    logmin=log10(min);
    bin_vuoti=0;

    while(count <= maxposition){
      if(log10(U->co[count]) < logmin+0.5*logdelta){
	while(log10(U->co[count]) < logmin+0.5*logdelta  && count<=maxposition){
	  count++;
	}

	bins->next=new_longpair();
	bins=bins->next;	
	bins->i=count-count1;
	/*			printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
				scanf("%c",&ch);
	*/
	bins->j=1;
	count1=count-1;
	(head->i)++;
	count1=count;
	count++;
      }else{
	bin_vuoti++;
	printf("\nWarning:: An empty bin was encountered\n");
      }			
      logmin+=logdelta;
    }
    if(bin_vuoti!=0){
      printf("\nThere was look for (%ld) empty bins\n",bin_vuoti);
    }
  }

  l=(DOUBLEBIN *)malloc(sizeof(DOUBLEBIN));
  if (!l) t_error("allocation failure in new_doublebin()");
  l->isdynamic=isDynamic;
  if(head->i < 2 ){
    t_error("Something wrong happened in binning");
  }else{
    l->index=new_longvector(head->i);
    (l->index)->nl=1;
    (l->index)->nh=head->i;
    bins=head->next;
    for(i=1;i<=head->i;i++){
      (l->index)->co[i]=bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }
    l->co=(double **)malloc((size_t) ((l->index->nh+NR_END)*sizeof(double *)));
    if (!l->co) t_error("allocation failure in new_doublebin()");
    l+=NR_END-1;

    l->co[1]=U->co+minposition-1;
    if(!(l->co[1])) t_error("Failure in splitting this vector");
    l->co[1]+=NR_END;
    l->co[1]-=1;
    bins=head->next;	
    for(i=2;i<=l->index->nh; i++){
      l->co[i]=l->co[i-1]+bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }

  }


  U->isdynamic=-1;
  U->nl=-1;
  U->nh=1;

  free(U);

  delete_longpair_list(head);

  return l;


}

/*---------------------------------------------------------------------------*/

double split2realvectors(DOUBLEVECTOR *U,DOUBLEVECTOR *T,DOUBLEBIN *l,DOUBLEBIN *m, long N,long num_max,FLOATVECTOR *novalue)


{

  double delta,min,max;
  long i,count,count1,minposition,maxposition,bin_vuoti;
  LONGPAIR *head,*bins;


  minposition=1;
  maxposition=U->nh;
  /*printf("%d\n",N);*/
  /*scanf("%c",&ch);*/
  if(N <=1){
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=1;
    count=2;
    while(count <= U->nh){
      while(U->co[count]==U->co[count-1] && count <=U->nh){
	count++;
      }
	
      bins->next=new_longpair();
      bins=bins->next;
      bins->i=count-count1;
      /*printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
	scanf("%c",&ch);*/

      bins->j=1;
      count1=count-1;
      head->i++;
      count1=count;
      count++; 
      if(head->i>num_max)t_error("Il numero di bin eccede il numero massimo che hai consentito");

    }

  }else if(N >1){
    if(novalue->co[1]==0){
      minposition=1;
      maxposition=U->nh;
    }else if(novalue->co[1]==-1){
      minposition=2;
      /*printf("*****%d\n",novalue);*/
      /*scanf("%hd",&ch);*/
      max=U->co[U->nh];
      while(U->co[minposition]<=novalue->co[2]){
	/*printf("%d\n",minposition);
	  scanf("%c",&ch);*/
	minposition++;
      }
      min=U->co[minposition];
      maxposition=U->nh;
    }else{
      maxposition=U->nh-1;
      while(U->co[maxposition]>=novalue->co[2]){
	/*printf("%d\n",maxposition);
	  scanf("%c",&ch);*/
	maxposition--;
      }	    
      max=U->co[maxposition];
      minposition=1;
      min=U->co[1];
    }

    printf("THE MINIMUM VALUE IS %f\n",min);
    printf("THE MAXIMUM VALUE IS %f\n",max);
    delta=(max-min)/(N-1);
    printf("DELTA IS%f\n",delta);
    /*  scanf("%c",&ch); */
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=minposition;
    count=minposition+1;
    bin_vuoti=0;
    while(count <= maxposition){
      if(U->co[count] < min +0.5*delta){
	while(U->co[count] < min+0.5*delta  && count<=maxposition){
	  count++;
	}
	bins->next=new_longpair();
	bins=bins->next;	
	bins->i=count-count1;
	/*printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
	  scanf("%c",&ch);*/

	bins->j=1;
	count1=count-1;
	(head->i)++;
	count1=count;
	count++; 

      }else{
	bin_vuoti++;
	printf("\nWarning:: An empty bin was encountered\n");
      }
      min+=delta;
    }
    if(bin_vuoti!=0){
      printf("\nThere was look for (%ld) empty bins\n",bin_vuoti);
    }

  }
  if(head->i < 2 ){
    t_error("Something wrong happened in binning");
  }else{
    m->index=new_longvector(head->i);
    l->index=new_longvector(head->i);
    (l->index)->nl=1;	
    (l->index)->nh=head->i;
    (m->index)->nl=1;	
    (m->index)->nh=head->i;

    bins=head->next;
    for(i=1;i<=head->i;i++){
      (l->index)->co[i]=bins->i;
      (m->index)->co[i]=bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }
    l->co=(double **)malloc((size_t) ((l->index->nh+NR_END)*sizeof(double *)));
    m->co=(double **)malloc((size_t) ((m->index->nh+NR_END)*sizeof(double *)));

    if (!l->co || !m->co) t_error("allocation failure in new_doublebin()");
    l+=NR_END-1;
    l->co[1]=U->co+minposition-1;
    m+=NR_END-1;
    m->co[1]=T->co+minposition-1;
		
    if(!(l->co[1]) || !(m->co[1])) t_error("Failure in splitting this vector");
    l->co[1]+=NR_END;
    l->co[1]-=1;
    m->co[1]+=NR_END;
    m->co[1]-=1;
		
    bins=head->next;	
    for(i=2;i<=l->index->nh; i++){
      l->co[i]=l->co[i-1]+bins->i;
      m->co[i]=m->co[i-1]+bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }

  }
  U->isdynamic=-1;
  U->nl=-1;
  U->nh=1;

  free(U);
  T->isdynamic=-1;
  T->nl=-1;
  T->nh=1;

  free(T);

  delete_longpair_list(head);

  if(N < 2) delta=0;
  return delta;

}



/*---------------------------------------------------------------------------*/

double esponentialsplit2realvectors(DOUBLEVECTOR *U,DOUBLEVECTOR *W, DOUBLEBIN* l,DOUBLEBIN* m,long N,long num_max,double base,FLOATVECTOR *novalue)


{

  double min,max,logdelta,logbase,logmin;
  long i,j,count,count1,minposition,maxposition,bin_vuoti;
  LONGPAIR *head,*bins;
  short mode;
  SHORTVECTOR *segna;
  minposition=1;
  maxposition=U->nh;
  logbase=log10(base);

  if(N <=1){
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=1;
    count=2;

    while(count <= U->nh){
      while(U->co[count]==U->co[count-1] && count <=U->nh){
	count++;
      }
	
      bins->next=new_longpair();
      bins=bins->next;
      bins->i=count-count1;
      /*printf("%f -> %d-%d=%d\n",U->co[count-1],count,count1,bins->i);
	scanf("%c",&ch);*/

      bins->j=1;
      count1=count-1;
      head->i++;
      count1=count;
      count++; 
      if(head->i>num_max)t_error("Il numero di bin eccede il numero massimo che hai consentito");
    }
  }else if(N >1){
    if(novalue->co[1]==0){
      minposition=1;
      if(U->co[minposition] <=0){
	while(U->co[minposition]<=0){
	  minposition++;
	  if( minposition > U->nh) t_error("No valid data for esponential binning");
	}	    
      }
      maxposition=U->nh;
    }else if(novalue->co[1]==-1){
      minposition=2;
      /* printf("*****%d\n",novalue->co[1]); */
      max=U->co[U->nh];
      while(U->co[minposition]<=novalue->co[2]){
	minposition++;
      }
      min=U->co[minposition];
      maxposition=U->nh;
    }else{
      maxposition=U->nh-1;
      while(U->co[maxposition]>=novalue->co[2]){
	/*			printf("%d\n",maxposition);
				scanf("%c",&ch);	
	*/			maxposition--;
      }	    
      max=U->co[maxposition];
      minposition=1;
      min=U->co[1];
    }
    printf("THE MINIMUM VALUE IS \t%f\n",min);
    printf("THE MAXIMUM VALUE IS \t%f\n",max);
    logdelta=(log10(max)-log10(min))/(logbase*(N-1));
    printf("LOG_10  DELTA IS \t%f\t",logdelta);
    /*scanf("%c",&ch); */
    head=new_longpair();
    head->i=0;
    bins=head;
    count1=minposition;
    count=minposition+1;
    logmin=log10(min);
    segna=new_shortvector(N);
    j=0;
    bin_vuoti=0;
    while(count <= maxposition){
      mode=1; 
      j++;
      segna->co[j]=0;
      if(log10(U->co[count])<=logmin+0.5*logdelta){
	while(log10(U->co[count]) <= logmin+0.5*logdelta  && count<=maxposition){
	  count++;
	}
	bins->next=new_longpair();
	bins=bins->next;	
	bins->i=count-count1;
	/*printf("%d %f -> %d-%d=%d\n",j,U->co[count-1],count,count1,bins->i);
	  bins->j=1;*/
	count1=count-1;
	(head->i)++;
	count1=count;
	count++;
      }else{
	bin_vuoti++;
	printf("\nWarning:: An empty bin was encountered\n");
      }
      logmin+=logdelta;
    }
    if(bin_vuoti!=0){
      printf("\nThere was look for (%ld) empty bins\n",bin_vuoti);
    }
  }
  if(head->i < 2 ){
    t_error("Something wrong happened in binning");
  }else{
    m->index=new_longvector(head->i);
    l->index=new_longvector(head->i);
    (l->index)->nl=1;	
    (l->index)->nh=head->i;
    (m->index)->nl=1;	
    (m->index)->nh=head->i;
    bins=head->next;
    for(i=1;i<=head->i;i++){
      (l->index)->co[i]=bins->i;
      (m->index)->co[i]=bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }
    l->co=(double **)malloc((size_t) ((l->index->nh+NR_END)*sizeof(double *)));
    m->co=(double **)malloc((size_t) ((m->index->nh+NR_END)*sizeof(double *)));

    if (!l->co || !m->co) t_error("allocation failure in new_doublebin()");
    l+=NR_END-1;
    l->co[1]=U->co+minposition-1;
    m+=NR_END-1;
    m->co[1]=W->co+minposition-1;
		
    if(!(l->co[1]) || !(m->co[1])) t_error("Failure in splitting this vector");
    l->co[1]+=NR_END;
    l->co[1]-=1;
    m->co[1]+=NR_END;
    m->co[1]-=1;
		
    bins=head->next;	
    for(i=2;i<=l->index->nh; i++){
      l->co[i]=l->co[i-1]+bins->i;
      m->co[i]=m->co[i-1]+bins->i;
      bins=bins->next;
      if(bins==NULL) break;
    }

  }

  U->isdynamic=-1;
  U->nl=-1;
  U->nh=1;

  free(U);
  W->isdynamic=-1;
  W->nl=-1;
  W->nh=1;

  free(W);

  delete_longpair_list(head);

  return logdelta;

}




DOUBLEMATRIX *shrink_doublematrix(DOUBLEMATRIX *data,FLOATVECTOR *Q)

{

  long i,count=0;


  DOUBLEMATRIX *newm;

  if(Q->co[1]<0){

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] > Q->co[2]){
	count++;
      }
    }
  }else if(Q->co[1]>0){

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] < Q->co[2]){
	count++;
      }
    }

  } else{

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] != Q->co[2]){
	count++;
      }
    }
	
	
  }


  newm=new_doublematrix(count,2);

  count=1;

  if(Q->co[1]<0){

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] > Q->co[2]){
	newm->co[count][1]=data->co[i][1];
	newm->co[count][2]=data->co[i][2];
	       
	count++;
      }
    }
  }else if(Q->co[1]>0){

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] < Q->co[2]){
	newm->co[count][1]=data->co[i][1];
	newm->co[count][2]=data->co[i][2];       
	count++;
      }
    }

  } else{

    for(i=1;i<=data->nrh;i++){

      if(data->co[i][2] != Q->co[2]){
	newm->co[count][1]=data->co[i][1];
	newm->co[count][2]=data->co[i][2];
	       
	count++;
      }
    }
	
	
  }


  return newm;

}

/* interpolating_function */ 

DOUBLEMATRIX *interpolating_function(DOUBLEMATRIX *cleandata)

{
  long i;
  DOUBLEMATRIX * W;
	
  W=new_doublematrix(cleandata->nrh,2);

  for(i=1;i<=cleandata->nrh-1;i++){
    W->co[i][1]=(cleandata->co[i+1][2]-cleandata->co[i][2])/(cleandata->co[i+1][1]-cleandata->co[i][1]);
    W->co[i][2]=cleandata->co[i][2]-W->co[i][1]*cleandata->co[i][1];
  }
	
  return W;
}

/* interpolate */

/*

double interpolate(double x,DOUBLEMATRIX *cleandata,DOUBLEMATRIX *W)

{
long count;
double y=0;    

count=W->co[W->nrh][1];

if(x< cleandata->co[1][1] || x>cleandata->co[cleandata->nrh][1] ){
t_error("Out of bound");
} else if ( x < cleandata->co[count][1]){
while(x < cleandata->co[count][1]){
count--;
}
W->co[W->nrh][1]=count;	
y=W->co[count][1]*x+W->co[count][2];		  	
}else{
while(x > cleandata->co[count][1]){
count++;
}
count--;	
W->co[W->nrh][1]=count;
y=W->co[count][1]*x+W->co[count][2];
}

return y;
}

*/

/*--------------------------------------------------------------------------*/

double interpolate(double x,DOUBLEMATRIX *cleandata,DOUBLEMATRIX *W)

{
  long count;
  double y=0;    

  count=W->nrh;

  if(x< cleandata->co[1][1] || x>cleandata->co[cleandata->nrh][1] ){
    t_error("Out of bound");
  } else if ( x <= cleandata->co[count][1]){
    while(x <= cleandata->co[count][1]){
      count--;
      if(count<=0){ 
	count++;
	break;
      }
    }
	
    y=W->co[count][1]*x+W->co[count][2];
  }  	

  return y;
}

/*--------------------------------------------------------------------------*/


/* interpolate_floatmatrix	*/

void interpolate_floatmatrix(FLOATMATRIX *matrice, float dt, float istante, FLOATVECTOR *vettore)

     /**
	Nuova subroutine che interpola da una matrice di dati float il valore all'istante scelto. 
	Inputs: matrice: 	matrice di float con i dati
	dt:			dt dati [s]
	istante:	istante in cui si vuole interpolare
	Outputs:vettore:	vettore di float con i dati per quell'istante
     */
{
  long i1,j,time;
  double resto;

  time=floor(istante/dt);

  if((istante/dt-time)==0){
    i1=time;
  }else{
    i1=time+1;
  }

  if(i1<1) t_error("tempo troppo basso");
  if(i1>=matrice->nrh) t_error("tempo troppo alto");
  resto=(istante-(i1-1)*dt)/dt;
  for(j=1;j<=matrice->nch;j++){
    vettore->co[j]=matrice->co[i1][j]+(matrice->co[i1+1][j]-matrice->co[i1][j])*resto;
  }
}

/*--------------------------------------------------------------------------*/

/* interpolate_doublematrix	*/

void interpolate_doublematrix(DOUBLEMATRIX *matrice, float dt, float istante, DOUBLEVECTOR *vettore)

     /**
	Nuova subroutine che interpola da una matrice di dati double il valore all'istante scelto. 
	Inputs: matrice: 	matrice di double con i dati
	dt:			dt dati [s]
	istante:	istante in cui si vuole interpolare
	Outputs:vettore:	vettore di double con i dati per quell'istante
     */

{
  long i1,j,time;
  double resto;

  time=floor(istante/dt);

  if((istante/dt-time)==0){
    i1=time;
  }else{
    i1=time+1;
  }

  if(i1<1){
    printf("istant=%f\n",istante); 
    t_error("tempo troppo basso");
  }
  if(i1>=matrice->nrh){ 
    printf("istant=%f,imax=%f\n",istante,(matrice->nrh-1)*dt);
    t_error("tempo troppo alto");
  }
  resto=(istante-(i1-1)*dt)/dt;
  for(j=1;j<=matrice->nch;j++){
    vettore->co[j]=matrice->co[i1][j]+(matrice->co[i1+1][j]-matrice->co[i1][j])*resto;
  }
}

/*--------------------------------------------------------------------------*/
double mean_doublematrix_column(DOUBLEMATRIX* net,long column)
{

  double mean;
  long i;

  mean=0;

  for(i=1;i<=net->nrh;i++){
    mean+=net->co[i][column];
  }

  return mean/=net->nrh;

}
/*--------------------------------------------------------------------------*/

double variance_doublematrix_column(DOUBLEMATRIX* net,long column,double mean)

{

  double variance;
  long i;

  variance=0;

  for(i=1;i<=net->nrh;i++){
    variance+=(net->co[i][column]-mean)*(net->co[i][column]-mean);
  }

  return variance/=net->nrh;

}

/*--------------------------------------------------------------------------*/

double approximate_2_multiple(double number,double div)

{

  return number-fabs(fmod(number,div));

}

/*-------- ricampiona ------------------------------------------------------------------*/

DOUBLEMATRIX *ricampiona(DOUBLEMATRIX *cleandata, float ti2, float dt2, long n2, float dt1)

{
  long i,j,k,nint;
  DOUBLEMATRIX *result, *coeff, *xy, *interval, *trapezio;
  long count;
  float perc;

  FLOATVECTOR *unused;
  unused=new_floatvector(1);
  initialize_floatvector(unused,0);



  /* matrix with x e y data */
  xy=new_doublematrix(cleandata->nrh,2);
  initialize_doublematrix(xy,0);

  /* assign xy x data */
  for(i=1;i<=cleandata->nrh;i++){
    xy->co[i][1]=cleandata->co[i][1];
  }
  /*doublematrix_control(xy,"xy.txt","xy",PRINT); */

  /* matrix with results */	
  result=new_doublematrix(n2,cleandata->nch);
  initialize_doublematrix(result,0);

  /* doublematrix_dem(result,unused, unused, "result.txt","doublematrix of result",NOPRINT);
     doublematrix_dem(cleandata,unused, unused, "result.txt","doublematrix of cleandata",NOPRINT); */
  printf("ti2=%f,dt2=%f,n2=%ld,dt1=%f \n",ti2,dt2,n2,dt1);

  /* assign results x- data */
  for(i=1;i<=n2;i++){
    /* result->co[i][1]=ti2+i*dt2-dt1; */
    result->co[i][1]=ti2+dt2/2.+(i-1)*dt2;
  } 
  /*doublematrix_control(result,"result.txt","result",PRINT);

  //matrix with intervals */
  interval=new_doublematrix(n2+1,2);
  initialize_doublematrix(interval,0);

  /* assign intervals x data: points where interpolate data*/
  /* interval->co[1][1]=ti2; */
  for(i=1;i<=n2+1;i++){
    /* interval->co[i][1]=ti2+(i-1)*dt2-dt1; */
    interval->co[i][1]=ti2+(i-1)*dt2;
	
  }
  /* doublematrix_control(interval,"prova.interval.txt","interval",PRINT); */

  /* matrix with point for each interval */
  trapezio=new_doublematrix(cleandata->nrh/n2+10,2);
  initialize_doublematrix(trapezio,0);

  /*doublematrix_control(trapezio,"prova.trapezio.txt","trapezio",PRINT);*/

  perc=0;
  count=0;

  /* make calculations for each cols of cleandata */
  for(j=2;j<=cleandata->nch;j++){

    /* assign xy y data */
    for(i=1;i<=cleandata->nrh;i++){
      xy->co[i][2]=cleandata->co[i][j];
    }

    /* creates a matrix with regession par */
    coeff=interpolating_function(xy);
    /*doublematrix_control(coeff,"prova.coeff.txt","coeff",PRINT);
      doublematrix_control(xy,"prova.xy.txt","xy",PRINT); */
	
    /* works for each interval */
    for(i=1;i<=n2;i++){
	
      /* writes controls */
		
      count++;
      if(count==100){
	perc=(i+(j-1)*n2)/((float)(cleandata->nch*n2))*100.;
	count=0;
	printf("computation %4.1f perc\n",perc);
      }

      /* count number of xy in the interval and assign x-values to trapezio */
      nint=1;
      trapezio->co[1][1]=interval->co[i][1];
      for(k=1;k<=xy->nrh;k++){
	if((xy->co[k][1]>interval->co[i][1])&&(xy->co[k][1]<interval->co[i+1][1])){
	  nint+=1;
	  trapezio->co[nint][1]=xy->co[k][1];
	}
      }
      nint+=1;
      trapezio->co[nint][1]=interval->co[i+1][1];
		
      /* assign y-values to trapezio */
      for(k=1;k<=nint;k++){
	trapezio->co[k][2]=interpolate(trapezio->co[k][1],xy,coeff);
      }
      /* correct bag of interpolate */
      if(trapezio->co[1][1]==xy->co[1][1]){
	trapezio->co[1][2]=xy->co[1][2];
      }
		
      /* calculate mean value to trapezio */
      result->co[i][j]=mean_function(trapezio, nint);
		
      /*printf("mean[%d][%d]=%f",i,j,result->co[i][j]);
	doublematrix_control(trapezio,"prova.trapezio.txt","trapezio",PRINT); */
    }
  }
	
  return result;
}


/*-------- quickinterpolate ------------------------------------------------------------------*/

DOUBLEMATRIX *quickinterpolate(DOUBLEMATRIX *cleandata, short nint)

     /* input: cleandata: double matrix with original data
	nint: number of intervals to aggegate
	return: double matrix with the interpolated data with (cleandata->nrh)/nint rows */

{

  long i,j,h,k,row2,row1;
  DOUBLEMATRIX *result;

  row1=cleandata->nrh;
  row2=row1/nint;

  result=new_doublematrix(row2,cleandata->nch);
  initialize_doublematrix(result,0);

  for(j=1;j<=result->nch;j++){
    for(i=1;i<=result->nrh;i++){
      for(k=1;k<=nint;k++){
	h=(i-1)*nint+k;
	result->co[i][j]+=cleandata->co[h][j];
      }
      result->co[i][j]/=(float)nint;		
    }
  }

  return result;
}


/*-------- mean_function ------------------------------------------------------------------*/

double mean_function(DOUBLEMATRIX *data, long n)

{

  double mean;
  long i;

  mean=0;

  for(i=1;i<n;i++){
    mean+=(data->co[i+1][2]+data->co[i][2])*
      (data->co[i+1][1]-data->co[i][1]);
  }


  return mean/=(2.*(data->co[n][1]-data->co[1][1]));

}

/*-------- cleandata_matrix ------------------------------------------------------------------*/


DOUBLEMATRIX *cleandata_matrix(DOUBLEMATRIX *data, FLOATVECTOR *V, SHORTMATRIX *control)

{
  DOUBLEMATRIX *cleandata;
  double nv, x_i, x_f;
  long i,j, n_i, n_f;

  nv=V->co[2];
  cleandata=new_doublematrix(data->nrh,data->nch);
  initialize_shortmatrix(control,0);
  copy_doublematrix(data,cleandata);

  for(j=1;j<=data->nch;j++){
    /* if nv in first row...*/
    i=1;
    n_i=0;
    n_f=0;
    x_i=0;
    x_f=0;
    while(i<data->nrh){
      while(data->co[i][j]!=nv && i<data->nrh) i++;
      x_i=data->co[i-1][j];
      n_i=i-1;
      if(i==data->nrh) n_i=0;
      if(n_i!=0){
	while(data->co[i][j]==nv && i<data->nrh) i++;
	x_f=data->co[i][j];
	n_f=i;
	/* call fill_data in datamanipulation.c */
	fill_data(cleandata,control,j,n_i,n_f,x_i,x_f);
			
			
      }
    }

    /* if nv in last row...*/
  }

  return cleandata;
}


/*-------- fill_data ------------------------------------------------------------------*/

void fill_data(DOUBLEMATRIX *cleandata,SHORTMATRIX *control,long j,long n_i,long n_f,double x_i,double x_f)

{

  double m;
  long i,n;

  n=n_f-n_i;
  m=(x_f-x_i)/n;

  for(i=1;i<n;i++){
    cleandata->co[n_i+i][j]=x_i+m*i;
    control->co[n_i+i][j]=1;
  }

}

/*-------- aggregate ------------------------------------------------------------------*/

DOUBLEMATRIX *aggregate(DOUBLEMATRIX *data, long col, float nv)

     /* calculates the mean aggregating all data wich have the same value in the column col
	input: data: double matrix with original data
	long col: the column with, 
	float nv: novalue
	return: double matrix aggregated data as mean values; 
	in the last column you have the number of aggregated elements */
 
{
  long i,j,nrows,count,num;
  DOUBLEMATRIX *result;
  FLOATVECTOR *U;
  U=new_floatvector(1);
  doublematrix_dem(data,U, U,"data.tmp", "data",NOPRINT);

  /* count rows */ 
  nrows=1;
  for(i=2;i<=data->nrh;i++){
    if(data->co[i-1][col]!=data->co[i][col]) nrows++;
  }
  result=new_doublematrix(nrows,data->nch+1);

  /* first row */
  for(j=1;j<=data->nch;j++){
    if(data->co[1][j]!=nv){
      result->co[1][j]=data->co[1][j];
    }else{
      result->co[1][j]=nv;
    }
  }
  result->co[1][data->nch+1]=1;
  /* makes mean */ 
  count=1;
  num=1;
  for(i=2;i<=data->nrh;i++){
    if(data->co[i-1][col]!=data->co[i][col]) {
      for(j=1;j<=data->nch;j++){
	if(result->co[count][j]!=nv){
	  result->co[count][j]/=num;
	}
	if(data->co[i][j]!=nv){
	  result->co[count+1][j]=data->co[i][j];
	}else{
	  result->co[count+1][j]=nv;
	}
      }
      result->co[count][data->nch+1]=num;
      count++;
      num=1;
		
    }else{
      for(j=1;j<=data->nch;j++){
	result->co[count][j]+=data->co[i][j];
      }
      num++;
    }
  }
  /* last row */
  if(result->co[count][j]!=nv){
    for(j=1;j<=data->nch;j++){
      result->co[count][j]/=num;
    }
  }
  result->co[count][data->nch+1]=num;

  doublematrix_dem(result,U, U,"result.tmp", "result",NOPRINT);
  return result;
}

