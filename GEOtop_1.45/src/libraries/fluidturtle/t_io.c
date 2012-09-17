//#include "turtle.h"
//#include "t_utilities.h"
//#include "tensor3D.h"
#include "t_io.h"
char *WORKING_DIRECTORY='\0';

/*WORKING_POSITION=SEEK_SET;*/

FILE * EXTERNAL_FILE;

LONGVECTOR *EXTERNAL_P;

char  EXTERNAL_FILE_NAME[256],OLD_NAME[256];

long int EXTERNAL_FILE_POSITION=SEEK_SET;

HEADER EXTERNAL_HEADER;



short OPENYES=0;

long IO_FILES_COUNTER=0;
long IO_STRINGS_COUNTER=0;
long  WORKING_POSITION=0;
long IO_PARMS_COUNTER=0;


t_keywords T_KEYWORDS={{"2","ascii","binary"},

							  {"7","char","short","int","long","float","double","string"},

							  {"5","array","vector","matrix","list","tensor"},

							  {"2","->","<-"},

							  {"2"," ",","},

							  {"2","{","}"}};




/**-----------------------------------------------------------------------*/
FILE *t_fopen(const char *name,const char *mode)
/* Safe file open */
{

//char newname[256];
FILE *fp=NULL;


/*if(strcmp(mode,"w")==0 || strcmp(mode,"wb")==0){
	if((fp=fopen(name,"r"))!=NULL ){
		// The file already exist
		//printf("\nWarning::Overwriting the file %s \n",name);
		strcpy(newname,name);
		strcat(newname,".old");
		t_fclose(fp);
		//rename(name,newname);
	}
}*/
    
if((fp=fopen(name,mode))==NULL){
	printf("%s",name);    
	t_error(" The specified file could not be opened ");

	return NULL;
	
}else{

	return fp;

}

}



/**-----------------------------------------------------------------------*/
FILE * t_fclose(FILE * stream)
/* Safe file close */
{


if(stream==NULL){
	printf(" An attemp was made to close an already closed file ");
}else{		
	 fclose(stream);
}

	 return NULL;

}




/**-----------------------------------------------------------------------*/
long read_shortvector_elements(FILE *input, SHORTVECTOR *v,char *mode)
/** 
It reads a vector of short  numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{


long count=0,i,tmp=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=v->nl;i<=v->nh;i++){
			tmp=fscanf(input,"%hd",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((short *)&(v->co[v->nl]),sizeof(short),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_shortvector_elements(FILE * output, SHORTVECTOR *v, long maxcols)
/** Write a vector of short int to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%d ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 		
 			}

	}
	putc('\n',output);
}

return OK;

}

/**-----------------------------------------------------------------------*/

long binarywrite_shortvector_elements(FILE * output, SHORTVECTOR *v)
/** Write a vector of short int to a file in binary mode */

{


long count=0,tmp=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((short *)&(v->co[v->nl]),sizeof(short),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_shortvector_elements(SHORTVECTOR *v,long maxcols)
/* Write a vector of short int to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%hd ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
putchar('\n');
}



}

/**-----------------------------------------------------------------------*/
long read_intvector_elements(FILE *input, INTVECTOR *v,char *mode)
/** 
It reads a vector of short  numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{

size_t i,tmp=0;
long count=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=v->nl;i<=v->nh;i++){
			tmp=fscanf(input,"%d",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((int *)&(v->co[v->nl]),sizeof(int),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_intvector_elements(FILE * output, INTVECTOR *v, long maxcols)
/** Write a vector of  int to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%d ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp ==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 		
 			}

	}
	putc('\n',output);
}

	return OK;



}


/**-----------------------------------------------------------------------*/

long binarywrite_intvector_elements(FILE * output, INTVECTOR *v)
/** Write a vector of  int to a file in binary mode */

{

long tmp=0;
long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((int *)&(v->co[v->nl]),sizeof(int),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_intvector_elements(INTVECTOR *v,long maxcols)
/* Write a vector of  int to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%d ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
	
}

putchar('\n');

}

/**-----------------------------------------------------------------------*/
long read_longvector_elements(FILE *input, LONGVECTOR *v,char *mode)
/** 
It reads a vector of long  numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{

size_t i,tmp=0;
long count=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=v->nl;i<=v->nh;i++){
			tmp=fscanf(input,"%ld",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((long *)&(v->co[v->nl]),sizeof(long),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_longvector_elements(FILE * output, LONGVECTOR *v, long maxcols)
/** Write a vector of  int to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%ld ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 		
 			}

	}
	putc('\n',output);
}

return OK;

}



/**-----------------------------------------------------------------------*/

long binarywrite_longvector_elements(FILE * output, LONGVECTOR *v)
/** Write a vector of  long to a file in binary mode */

{

long tmp=0;
long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((long *)&(v->co[v->nl]),sizeof(long),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_longvector_elements(LONGVECTOR *v,long maxcols)
/* Write a vector of  long to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%ld ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
	
}

putchar('\n');

}

/**-----------------------------------------------------------------------*/
long read_floatvector_elements(FILE *input, FLOATVECTOR *v,char *mode)
/** 
It reads a vector of  float numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{

size_t i,tmp=0;
long count=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=v->nl;i<=v->nh;i++){
			tmp=fscanf(input,"%f",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((float *)&(v->co[v->nl]),sizeof(float),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_floatvector_elements(FILE * output, FLOATVECTOR *v, long maxcols)
/** Write a vector of  float to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%f ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF; 		
 			}

	}
	putc('\n',output);
}

return OK;

}



/**-----------------------------------------------------------------------*/

long binarywrite_floatvector_elements(FILE * output, FLOATVECTOR *v)
/** Write a vector of  float to a file in binary mode */

{

long tmp=0;
long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((float *)&(v->co[v->nl]),sizeof(float),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_floatvector_elements(FLOATVECTOR *v,long maxcols)
/* Write a vector of float to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%f ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
	
}

putchar('\n');

}

/**-----------------------------------------------------------------------*/
long read_doublevector_elements(FILE *input, DOUBLEVECTOR *v,char *mode)
/** 
It reads a vector of  float numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{

size_t i,tmp=0;
long count=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=v->nl;i<=v->nh;i++){
			tmp=fscanf(input,"%lf",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((double *)&(v->co[v->nl]),sizeof(double),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_doublevector_elements(FILE * output, DOUBLEVECTOR *v, long maxcols)
/** Write a vector of  double to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%f ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 		
 			}

	}
	putc('\n',output);
}

return OK;

}



/**-----------------------------------------------------------------------*/

long binarywrite_doublevector_elements(FILE * output, DOUBLEVECTOR *v)
/** Write a vector of  double to a file in binary mode */

{

long tmp=0;
long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((double *)&(v->co[v->nl]),sizeof(double),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_doublevector_elements(DOUBLEVECTOR *v,long maxcols)
/* Write a vector of double to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%f ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
	
}

putchar('\n');

}

/**-----------------------------------------------------------------------*/
long read_charvector_elements(FILE *input, CHARVECTOR *v,char *mode)
/** 
It reads a vector of  char numbers and returns the number of
elements read. Storage mode can be either ascii or binary.
*/

{

size_t i,tmp=0;
long count=0;
const char ascii[2]="a",binary[2]="b";

if(input==NULL){
	t_error("The input file was not opened properly or is not allocated");
}else if (v==NULL || v->co==NULL || (v->isdynamic)!=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh){
	t_error("The vector has no proper dimensions");
} else if(strcmp(mode,ascii)==0){
	for(i=v->nl;i<=v->nh;i++){
	        v->co[i]='\0';
			tmp=fscanf(input,"%c",&(v->co[i]));
			if(tmp!=EOF){
				count+=tmp;
			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf(" after position %ld\n",count);
				return -count;
			}
	}
} else if(strcmp(mode,binary)){
	count=fread((char *)&(v->co[v->nl]),sizeof(char),v->nh-v->nl+1,input);
}else{
	t_error("Error in reading mode::Mode not supported");
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_charvector_elements(FILE * output, CHARVECTOR *v, long maxcols)
/** Write a vector of  char to a file */

{

	long i,tmp=0;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			tmp=fprintf(output,"%c ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putc('\n',output);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 		
 			}

	}
	putc('\n',output);
}

return OK;

}



/**-----------------------------------------------------------------------*/

long binarywrite_charvector_elements(FILE * output, CHARVECTOR *v)
/** Write a vector of  char to a file in binary mode */

{

long tmp=0;
long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {
 	tmp=fwrite((char *)&(v->co[v->nl]),sizeof(char),v->nh-v->nl+1,output);
 	if(tmp!=EOF){
 		count+=tmp;
 	}else{
		printf("Error in stored data::Unespected End of file encountered");
		printf("after position %ld\n",count);
		return -count;
 		
 	}
}

if(count!=(v->nh-v->nl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}


/**-----------------------------------------------------------------------*/
void print_charvector_elements(CHARVECTOR *v,long maxcols)
/* Write a vector of char to the standard output */

{

long i;

putchar('\n');

if (v==NULL || v->co==NULL || v->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(v->nl > v->nh ){
	t_error("The vector has no proper dimensions");
} else {

	for(i=v->nl;i<=v->nh;i++){
			printf("%c ",v->co[i]);
		    if(i%maxcols==0 && i!=(v->nh)) putchar('\n');
	}
		
	
}

putchar('\n');

}

/**----------------------------------------------------------------------- 

Matrixes

-----------------------------------------------------------------------*/

long read_shortmatrix_elements(FILE *input,SHORTMATRIX *m,char *mode)

/** Read a matrix of short  numbers */

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j;



if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch){
	t_error("The matrix has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fscanf(input,"%hd",&(m->co[i][j]));
			if(tmp!=EOF){
				count+=tmp;
			}else {
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
				return -count;
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
  	
  	for(i=m->nrl;i<=m->nrh;i++){
  		tmp=fread((short *)&(m->co[i][m->ncl]),sizeof(short),m->nch-m->ncl+1,input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_shortmatrix_elements(FILE * output,SHORTMATRIX *m,long maxcols)
/* Write a matrix of floating point numbers to a file */

{

	long tmp=0,i,j;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fprintf(output,"%hd\t ",m->co[i][j]);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 			}

		    if(j%maxcols==0 && j!=m->nch) putc('\n',output);
		}
		putc('\n',output);
	}
	putchar('\n');
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_shortmatrix_elements(FILE * output, SHORTMATRIX *m)
/** Write a matrix of shorts  to a file */

{


long tmp,i,count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(m->ncl > m->nch || m->nrl > m->nrh){
	t_error("The matrix has no proper dimensions");
} else {
	for(i=m->nrl;i<=m->nrh;i++){
 		tmp=fwrite((short *)&(m->co[i][m->nrl]),sizeof(short),m->nch-m->ncl+1,output);
 		if(tmp!=EOF){
 			count+=tmp;
 		}else{
			printf("Error in stored data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
 		
 		}
 	}
}

if(count!=(m->nch-m->ncl+1)*(m->nrh-m->nrl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}

/**-----------------------------------------------------------------------*/
void print_shortmatrix_elements(SHORTMATRIX *m,long maxcols)
/* Write a matrix of  shorts to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%hd ",m->co[i][j]);
		    if(j%maxcols==0 && j!=m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');
}


}
/*-----------------------------------------------------------------------*/

long read_intmatrix_elements(FILE *input,INTMATRIX *m,char *mode)
/** 

Read a matrix of int numbers

*/

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j;



if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch){
	t_error("The matrix has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fscanf(input,"%d",&(m->co[i][j]));
			if(tmp!=EOF){
				count+=tmp;
			}else {
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
				return -count;
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
  	
  	for(i=m->nrl;i<=m->nrh;i++){
  		tmp=fread((int *)&(m->co[i][m->ncl]),sizeof(int),m->nch-m->ncl+1,input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_intmatrix_elements(FILE * output,INTMATRIX *m,long maxcols)
/* Write a matrix of int numbers to a file */

{

	long tmp=0,i,j;
	//long count=0;

if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fprintf(output,"%d ",m->co[i][j]);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 			}

		    if(j%maxcols==0 && j!=m->nch) putc('\n',output);
		}
		putc('\n',output);
	}
	putchar('\n');
}

return OK;
}



/**-----------------------------------------------------------------------*/
long binarywrite_intmatrix_elements(FILE * output, INTMATRIX *m)
/** Write a matrix of int  to a file */

{

long tmp=0;
long i,count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(m->ncl > m->nch || m->nrl > m->nrh){
	t_error("The matrix has no proper dimensions");
} else {
	for(i=m->nrl;i<=m->nrh;i++){
 		tmp=fwrite((int *)&(m->co[i][m->nrl]),sizeof(int),m->nch-m->ncl+1,output);
 		if(tmp!=EOF){
 			count+=tmp;
 		}else{
			printf("Error in stored data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
 		
 		}
 	}
}

if(count!=(m->nch-m->ncl+1)*(m->nrh-m->nrl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}

/**-----------------------------------------------------------------------*/
void print_intmatrix_elements(INTMATRIX *m,long maxcols)
/* Write a matrix of  double to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%d ",m->co[i][j]);
		    if(j%maxcols==0 && j!=m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');
}


}



/*-----------------------------------------------------------------------*/

long read_longmatrix_elements(FILE *input,LONGMATRIX *m,char *mode)
/** Read a matrix of long point numbers */

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j;

if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch){
	t_error("The matrix has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fscanf(input,"%ld",&(m->co[i][j]));
/*			fscanf(input,"%c",&ch);
			printf("^^^%d %c\n",m->co[i][j],ch);
			scanf("%c",&ch);
*/			if(tmp!=EOF){
				count+=tmp;
/*				printf("^@^%d\n",count); */
			}else {
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
				return -count;
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
  	
  	for(i=m->nrl;i<=m->nrh;i++){
  		tmp=fread((long *)&(m->co[i][m->ncl]),sizeof(long),m->nch-m->ncl+1,input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)){
/*    printf("^^^%d %d\n",count,(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)); */
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_longmatrix_elements(FILE * output,LONGMATRIX *m,long maxcols)
/* Write a matrix of long point numbers to a file */

{

	long tmp=0,i,j;
	//long count=0;


if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			/*printf("i(%d,%d)m(%d)\n",i,j,m->co[i][j]);*/
			tmp=fprintf(output,"%ld ",m->co[i][j]);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 			}

		    if(j%maxcols==0 && j!=m->nch) putc('\n',output);
		}
		putc('\n',output);
	}
}

return OK;
}



/**-----------------------------------------------------------------------*/
long binarywrite_longmatrix_elements(FILE * output, LONGMATRIX *m)
/** Write a matrix of long  to a file */

{

long tmp=0;
long i,count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(m->ncl > m->nch || m->nrl > m->nrh){
	t_error("The matrix has no proper dimensions");
} else {
	for(i=m->nrl;i<=m->nrh;i++){
 		tmp=fwrite((long *)&(m->co[i][m->nrl]),sizeof(long),m->nch-m->ncl+1,output);
 		if(tmp!=EOF){
 			count+=tmp;
 		}else{
			printf("Error in stored data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
 		
 		}
 	}
}

if(count!=(m->nch-m->ncl+1)*(m->nrh-m->nrl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}

/**-----------------------------------------------------------------------*/
void print_longmatrix_elements(LONGMATRIX *m,long maxcols)
/* Write a matrix of  double to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%ld ",m->co[i][j]);
		    if(j%maxcols==0 && j!=m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');
}


}






/*-----------------------------------------------------------------------*/

long read_floatmatrix_elements(FILE *input,FLOATMATRIX *m,char *mode)
/** 

Read a matrix of floating point numbers

*/

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j;



if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch){
	t_error("The matrix has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fscanf(input,"%f",&(m->co[i][j]));
			if(tmp!=EOF){
				count+=tmp;
			}else {
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
				return -count;
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
  	
  	for(i=m->nrl;i<=m->nrh;i++){
  		tmp=fread((float *)&(m->co[i][m->ncl]),sizeof(float),m->nch-m->ncl+1,input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_floatmatrix_elements(FILE * output,FLOATMATRIX *m,long maxcols)
/* Write a matrix of floating point numbers to a file */

{

	long tmp=0,i,j;
	//long count=0;

if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fprintf(output,"%f ",m->co[i][j]);
			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 			}

		    if(j%maxcols==0 && j!=m->nch) putc('\n',output);
		}
		putc('\n',output);
	}
	putchar('\n');
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_floatmatrix_elements(FILE * output, FLOATMATRIX *m)
/** Write a matrix of float  to a file */

{

long tmp=0;
long i,count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(m->ncl > m->nch || m->nrl > m->nrh){
	t_error("The matrix has no proper dimensions");
} else {
	for(i=m->nrl;i<=m->nrh;i++){
 		tmp=fwrite((float *)&(m->co[i][m->nrl]),sizeof(float),m->nch-m->ncl+1,output);
 		if(tmp!=EOF){
 			count+=tmp;
 		}else{
			printf("Error in stored data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
 		
 		}
 	}
}

if(count!=(m->nch-m->ncl+1)*(m->nrh-m->nrl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}

/**-----------------------------------------------------------------------*/
void print_floatmatrix_elements(FLOATMATRIX *m,long maxcols)
/* Write a matrix of  float to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%f ",m->co[i][j]);
		    if(j%maxcols==0 && j!=m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');
}


}


/*-----------------------------------------------------------------------*/

long read_doublematrix_elements(FILE *input,DOUBLEMATRIX *m,char *mode)
/** Read a matrix of double numbers */

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j;



if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch){
	t_error("The matrix has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fscanf(input,"%lf",&(m->co[i][j]));
			if(tmp!=EOF){
				count+=tmp;
			}else {
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
				return -count;
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
  	
  	for(i=m->nrl;i<=m->nrh;i++){
  		tmp=fread((double *)&(m->co[i][m->ncl]),sizeof(double),m->nch-m->ncl+1,input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered\n");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_doublematrix_elements(FILE * output,DOUBLEMATRIX *m,long maxcols)
/* Write a matrix of double numbers to a file */

{

	long tmp=0,i,j;
	//long count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			tmp=fprintf(output,"%g\t",m->co[i][j]);
			if(tmp==EOF){
 				printf("Error in storing data::Unespected End of file encountered\n");
				return EOF;
 			}

		    if(j%maxcols==0 && j!= m->ncl) putc('\n',output);
		}
		putc('\n',output);
	}
	/*putchar('\n');*/
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_doublematrix_elements(FILE * output, DOUBLEMATRIX *m)
/** Write a matrix of double  to a file */

{


long i,tmp=0,count=0;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The vector was not allocated properly");
}else if(m->ncl > m->nch || m->nrl > m->nrh){
	t_error("The matrix has no proper dimensions");
} else {
	for(i=m->nrl;i<=m->nrh;i++){
 		tmp=fwrite((double *)&(m->co[i][m->nrl]),sizeof(double),m->nch-m->ncl+1,output);
 		if(tmp!=EOF){
 			count+=tmp;
 		}else{
			printf("Error in stored data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
 		
 		}
 	}
}

if(count!=(m->nch-m->ncl+1)*(m->nrh-m->nrl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}


}

/**-----------------------------------------------------------------------*/
void print_doublematrix_elements(DOUBLEMATRIX *m,long maxcols)
/* Write a matrix of  double to the standard output */

{

long i,j;

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
			printf("%f ",m->co[i][j]);
		    if(j%maxcols==0 && j!= m->nch) putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');	
}


}


/**-----------------------------------------------------------------------*/
long read_intbin_elements(FILE *input,INTBIN *l,char *mode)
/* Read a list of int  numbers */


{

long i,j,tmp=0,count=0,chksum=0;
const char ascii[2]="a",binary[2]="b";


if(input==NULL){
	t_error("The input file was not opened properly");
} else if(l==NULL  || l->co==NULL || (l->isdynamic)!=1    || l->index==NULL ){
	t_error("The bin was not allocated properly");

} else if(strcmp(mode,ascii)==0){
	
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
	
			tmp=fscanf(input,"%d",&(l->co[i][j]));
 			if(tmp!=EOF){
 				count+=tmp;
 			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf("after position %ld\n",count);
				return -count;

		
			}

		}
	
	}
  

}else if(strcmp(mode,binary)==0){

  	for(i=1;i<=(l->index)->nh;i++){
  		chksum+=l->index->co[i];
  		tmp=fread((int *)&(l->co[i][1]),sizeof(int),l->index->co[i],input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}


} else {

	    t_error("Error in reading mode::Mode not supported");
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_intbin_elements(FILE * output, INTBIN *l,long maxcols)
/** Write a list of int to a file */

{


	long tmp=0,i,j;
	//long count=0;
	long chksum=0;



if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
			tmp=fprintf(output,"%d ",l->co[i][j]);
 			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered");
				return EOF;			
			}
		    if(j%maxcols==0 && j!=(l->index)->co[i]) putc('\n',output);
		}
		putc('\n',output);
	}
	
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_intbin_elements(FILE * output, INTBIN *l)
/** Write a list of int to a file */

{



long i,tmp=0,count=0,chksum=0;

/** Here is still missing a check  that output was open as "wb" */

if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		tmp=fwrite((int *)&(l->co[i][1]),sizeof(int),l->index->co[i],output);
 		if(tmp!=EOF){	
			count+=tmp;
		}else {
			printf("Error in storing data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
			
		}
	}
	
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}





/**-----------------------------------------------------------------------*/
void print_intbin_elements(INTBIN *l,long maxcols)
/* Write a list of  int to the standard output */

{



long i,j;



if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putchar('\n');
	for(i=1;i<=(l->index)->nh;i++){
		for(j=1;j<=(l->index)->co[i];j++){
		/* This is the specialized code for int */
			printf("%d ",l->co[i][j]);
			if(j%maxcols==0 && j!=(l->index)->co[i]) putchar('\n');
		}
		putchar('\n');
	}
	
}
	
}

/**-----------------------------------------------------------------------*/
long read_shortbin_elements(FILE *input,SHORTBIN *l,char *mode)
/* Read a list of int  numbers */


{

long i,j,tmp=0,count=0,chksum=0;
const char ascii[2]="a",binary[2]="b";


if(input==NULL){
	t_error("The input file was not opened properly");
} else if(l==NULL  || l->co==NULL || (l->isdynamic)!=1    || l->index==NULL ){
	t_error("The bin was not allocated properly");

} else if(strcmp(mode,ascii)==0){
	
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
	
			tmp=fscanf(input,"%hd",&(l->co[i][j]));
 			if(tmp!=EOF){
 				count+=tmp;
 			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf("after position %ld\n",count);
				return -count;

		
			}

		}
	
	}
  

}else if(strcmp(mode,binary)==0){

  	for(i=1;i<=(l->index)->nh;i++){
  		chksum+=l->index->co[i];
  		tmp=fread((short *)&(l->co[i][1]),sizeof(short),l->index->co[i],input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}


} else {

	    t_error("Error in reading mode::Mode not supported");
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_shortbin_elements(FILE * output, SHORTBIN *l,long maxcols)
/** Write a list of int to a file */

{



	long tmp=0,i,j;
	//long count=0;
	long chksum=0;
	


if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
			tmp=fprintf(output,"%hd ",l->co[i][j]);
 			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered");
				return EOF;			
			}
		    if(j%maxcols==0 && j!=(l->index)->co[i]) putc('\n',output);
		}
		putc('\n',output);
	}
	
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_shortbin_elements(FILE * output, SHORTBIN *l)
/** Write a list of short to a file */

{



long i,tmp=0,count=0,chksum=0;

/** Here is still missing a check  that output was open as "wb" */

if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		tmp=fwrite((short *)&(l->co[i][1]),sizeof(short),l->index->co[i],output);
 		if(tmp!=EOF){	
			count+=tmp;
		}else {
			printf("Error in storing data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
			
		}
	}
	
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}





/**-----------------------------------------------------------------------*/
void print_shortbin_elements(SHORTBIN *l,long maxcols)
/* Write a list of  int to the standard output */

{



long i,j;



if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putchar('\n');
	for(i=1;i<=(l->index)->nh;i++){
		for(j=1;j<=(l->index)->co[i];j++){
		/* This is the specialized code for int */
			printf("%d ",l->co[i][j]);
			if(j%maxcols==0 && j!=(l->index)->co[i]) putchar('\n');
		}
		putchar('\n');
	}
	
}
	
}

/**-----------------------------------------------------------------------*/
long read_longbin_elements(FILE *input,LONGBIN *l,char *mode)
/* Read a list of long  numbers */


{

long i,j,tmp=0,count=0,chksum=0;
const char ascii[2]="a",binary[2]="b";


if(input==NULL){
	t_error("The input file was not opened properly");
} else if(l==NULL  || l->co==NULL || (l->isdynamic)!=1    || l->index==NULL ){
	t_error("The bin was not allocated properly");

} else if(strcmp(mode,ascii)==0){
	
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
	
			tmp=fscanf(input,"%ld",&(l->co[i][j]));
 			if(tmp!=EOF){
 				count+=tmp;
 			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf("after position %ld\n",count);
				return -count;

		
			}

		}
	
	}
  

}else if(strcmp(mode,binary)==0){

  	for(i=1;i<=(l->index)->nh;i++){
  		chksum+=l->index->co[i];
  		tmp=fread((long *)&(l->co[i][1]),sizeof(long),l->index->co[i],input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}


} else {

	    t_error("Error in reading mode::Mode not supported");
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_longbin_elements(FILE * output, LONGBIN *l,long maxcols)
/** Write a list of int to a file */

{



	long tmp=0,i,j;
	//long count=0;
	long chksum=0;
	


if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
			tmp=fprintf(output,"%ld ",l->co[i][j]);
 			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered");
				return EOF;			
			}
		    if(j%maxcols==0 && j!=(l->index)->co[i]) putc('\n',output);
		}
		putc('\n',output);
	}
	
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_longbin_elements(FILE * output, LONGBIN *l)
/** Write a list of long to a file */

{



long i,tmp=0,count=0,chksum=0;

/** Here is still missing a check  that output was open as "wb" */

if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		tmp=fwrite((long *)&(l->co[i][1]),sizeof(long),l->index->co[i],output);
 		if(tmp!=EOF){	
			count+=tmp;
		}else {
			printf("Error in storing data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
			
		}
	}
	
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}





/**-----------------------------------------------------------------------*/
void print_longbin_elements(LONGBIN *l,long maxcols)
/* Write a list of  long to the standard output */

{



long i,j;



if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putchar('\n');
	for(i=1;i<=(l->index)->nh;i++){
		for(j=1;j<=(l->index)->co[i];j++){
		/* This is the specialized code for long */
			printf("%ld ",l->co[i][j]);
			if(j%maxcols==0 && j!=(l->index)->co[i]) putchar('\n');
		}
		putchar('\n');
	}
	
}
	
}

/**-----------------------------------------------------------------------*/
long read_doublebin_elements(FILE *input,DOUBLEBIN *l,char *mode)
/* Read a list of double  numbers */


{

long i,j,tmp=0,count=0,chksum=0;
const char ascii[2]="a",binary[2]="b";


if(input==NULL){
	t_error("The input file was not opened properly");
} else if(l==NULL  || l->co==NULL || (l->isdynamic)!=1    || l->index==NULL ){
	t_error("The bin was not allocated properly");

} else if(strcmp(mode,ascii)==0){
	
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
	
			tmp=fscanf(input,"%lf",&(l->co[i][j]));
 			if(tmp!=EOF){
 				count+=tmp;
 			}else{
				printf("Error in stored data::Unespected End of file encountered");
				printf("after position %ld\n",count);
				return -count;

		
			}

		}
	
	}
  

}else if(strcmp(mode,binary)==0){

  	for(i=1;i<=(l->index)->nh;i++){
  		chksum+=l->index->co[i];
  		tmp=fread((double *)&(l->co[i][1]),sizeof(double),l->index->co[i],input);
  		if(tmp!=EOF){
  			count+=tmp;
  		}else{
				printf("Error in stored data::Unexpected End of File encountered");
				printf(" after position %ld\n",count);
  				return -count;
  		}
	}


} else {

	    t_error("Error in reading mode::Mode not supported");
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}

/**-----------------------------------------------------------------------*/
long write_doublebin_elements(FILE * output, DOUBLEBIN *l,long maxcols)
/** Write a list of double to a file */

{



	long tmp=0,i,j;
	//long count=0;
	long chksum=0;
	


if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
			tmp=fprintf(output,"%lf ",l->co[i][j]);
 			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered");
				return EOF;			
			}
		    if(j%maxcols==0 && j!=(l->index)->co[i]) putc('\n',output);
		}
		putc('\n',output);
	}
	
}

return OK;

}



/**-----------------------------------------------------------------------*/
long binarywrite_doublebin_elements(FILE * output, DOUBLEBIN *l)
/** Write a list of double to a file */

{



long i,tmp=0,count=0,chksum=0;

/** Here is still missing a check  that output was open as "wb" */

if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		tmp=fwrite((double *)&(l->co[i][1]),sizeof(double),l->index->co[i],output);
 		if(tmp!=EOF){	
			count+=tmp;
		}else {
			printf("Error in storing data::Unespected End of file encountered");
			printf("after position %ld\n",count);
			return -count;
			
		}
	}
	
}

if(count!=chksum){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
} else {
	return count;
}

}





/**-----------------------------------------------------------------------*/
void print_doublebin_elements(DOUBLEBIN *l,long maxcols)
/* Write a list of  double to the standard output */

{



long i,j;



if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putchar('\n');
	for(i=1;i<=(l->index)->nh;i++){
		for(j=1;j<=(l->index)->co[i];j++){
		/* This is the specialized code for double */
			printf("%f ",l->co[i][j]);
			if(j%maxcols==0 && j!=(l->index)->co[i]) putchar('\n');
		}
		putchar('\n');
	}
	
}
	
}



/**-----------------------------------------------------------------------*/


STRINGBIN *read_plane_strings(FILE *inputfile,long length,long maxbuffersize)

{

//short ind=1;

long count=0, ll=0, df;

char *buffer=NULL;

long int buffer_index=0, buffer_size=0;

LONGVECTOR *v;

STRINGBIN *m;



v=new_longvector(length);

skip_whitespaces(inputfile);

buffer_size=BUFFERINCREMENT;

		buffer=(char *)malloc((buffer_size+1)*sizeof(char));

		if(buffer==NULL){

			t_error("Cannot allocate the buffer");

		}

/* The buffer is NULL terminated */
		
buffer[buffer_size]='\0';

do{
		if(buffer_index < buffer_size){

			buffer[buffer_index]=fgetc(inputfile);
			
			if(isspace(buffer[buffer_index])){
			
			    count++;
				
				/* dimensions of the index of string bins is increased by one
				in order to contain the termination character */
				
				v->co[count]=ll+1;
				
				ll=0;
								
				buffer_index++; 
				
				skip_whitespaces(inputfile);


			}else {

			 ++buffer_index;
			 ll++;


			}

		}else {
	
		
			if(buffer_size==maxbuffersize){

				printf("Warning::A very long string has exceeded the maximum  buffer size\n");

				count++;
				
				v->co[count]=ll+1;
				
				ll=0;
				
				df=v->nh-count;
				
				printf ("Warning::missing part of string %ld and other %ld strings\n",count,df);
				
				v->nh=count;
				
				do{
			
					buffer[buffer_size-1]=getc(inputfile);
			
				}while(buffer[buffer_size-1]!=EOF);
				
				
				buffer[buffer_size-1]=' ';
							
						
				break; 
			
				
		
		 } else if(buffer_size + BUFFERINCREMENT > maxbuffersize){

			buffer_size=maxbuffersize;

			buffer=(char *)realloc(buffer,(buffer_size+1)*sizeof(char));

			if(buffer==NULL){

				t_error("Cannot expand the buffer");

			}
			
			buffer[buffer_size]='\0';
			
		} else {
		
			buffer_size+=BUFFERINCREMENT;

			buffer=(char *)realloc(buffer,(buffer_size+1)*sizeof(char));

			if(buffer==NULL){

				t_error("Cannot expand the buffer");

			}
			
			buffer[buffer_size]='\0';

		}
	
	
	}


	} while(count < length && buffer[buffer_index-1]!=EOF);

	if(buffer[buffer_index-1]==EOF){

		t_error("Unespected  End of file");

	}


m=new_stringbin(v);

copy_buffer_into_stringbin(buffer,m);

free_longvector(v);

free(buffer);


return m;

}


/**-----------------------------------------------------------------------*/
long write_stringbin_elements(FILE * output, STRINGBIN *l,long maxcols)
/** Write a list  of strings to a file */

{



	long tmp=0,i,j;
	//long count=0;
	long chksum=0;
	

if(output==NULL){
	t_error("The output file was not opened properly");
} else if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
} else {
	putc('\n',output);
	for(i=1;i<=(l->index)->nh;i++){
		chksum+=l->index->co[i];
		for(j=1;j<=(l->index)->co[i];j++){
			tmp=fprintf(output,"%c",l->co[i][j]);
 			if(tmp==EOF){
				printf("Error in storing data::Unespected End of file encountered");
				return EOF;			
			}
		    if(j%maxcols==0 && j!=(l->index)->co[i]-1) putc('\n',output);

		}
		putc('\n',output);
	}
putc('\n',output);	
}

return OK;

}



/**-----------------------------------------------------------------------*/
void print_stringbin_elements(STRINGBIN *l,long maxcols)
/* Write a list of  strings to the standard output */

{



long i,j;



if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
        //printf("%d  %d  %d\n",l,l->co,l->isdynamic);
        print_longvector_elements(l->index,100);
	t_error("The bin was not allocated properly");
} else {
	putchar('\n');
	for(i=1;i<=(l->index)->nh;i++){
	/* The index below do not go to the last - \0 - element */
		for(j=1;j<(l->index)->co[i];j++){
		/* This is the specialized code for characters */
			printf("%c",l->co[i][j]);
			if(j%maxcols==0 && j!=(l->index)->co[i]-1) putchar('\n');
		}
		putchar('\n');
	}
	
}
	
}


/**-----------------------------------------------------------------------*/
int copy_buffer_into_stringbin(char *buffer,STRINGBIN *l)
/* Write a buffer of strings  into its stringbin. "Its" means that
the string bin is already correctly allocated and inizialized */

{



long i,j,count=0;


if(l==NULL || l->co==NULL || (l->isdynamic)!=1  || l->index==NULL){
	t_error("The bin was not allocated properly");
	
}else if(buffer==NULL){
	t_error("The buffer is NULL");
} else {
	for(i=1;i<=(l->index)->nh;i++){
		for(j=1;j<(l->index)->co[i];j++){
		    if(buffer[count]!='\0'){
		    	l->co[i][j]=buffer[count];
             	           count++;
                          } else {
             	           t_error("Unespected end of buffer::The buffer has been scanned completely ");
                          }
		}
		if(isspace(buffer[count])) {
			count++;
		}
                     l->co[i][(l->index)->co[i]]='\0';
	}
	
}

return OK;
	
}




/**-----------------------------------------------------------------------*/

long longvectorcmp(LONGVECTOR *A,LONGVECTOR *B)

{

long count=0,i;

if(A->isdynamic!=B->isdynamic) {

count=1;

}

if(A->nh!=B->nh) {

count +=10;

}


for(i=1;i<=A->nh; i++){

	if(A->co[i]!=B->co[i]){
	
		count+=100;		
	
	}

}

return count;


}

/**-----------------------------------------------------------------------*/

void skip_whitespaces(FILE *inputfile)

{

int whitespace;
char ctmp;

do{
	whitespace=fgetc(inputfile);
	if(whitespace==EOF) printf("\nWarning::Unespected End of File encountered\n");
}while(isspace(whitespace) || iscntrl(whitespace));

ctmp=ungetc(whitespace,inputfile);

}

/**-----------------------------------------------------------------------*/

void goto_EOF(FILE *inputfile)



{



int whitespace;
char ctmp;


	do{

		whitespace=fgetc(inputfile);
		
			if(whitespace==EOF) break;

		}while(isspace(whitespace));

	ctmp=ungetc(whitespace,inputfile);




}




/**-----------------------------------------------------------------------*/

int iscomment(char *buffer, FILE *inputfile)

{


long buffer_index=0, pos;


	 pos=ftell(inputfile);
	 skip_whitespaces(inputfile);
  	 while(buffer_index<3){

		buffer[buffer_index]=fgetc(inputfile);
        buffer_index++;

	}

	buffer[buffer_index]='\0'; /* buffer_index is 3 */


	if(strcmp(buffer,"/**")!=0) {
		 fseek(inputfile,pos,SEEK_SET);	
		 return buffer_index=0;

	} else {

		return buffer_index;

	}

}

/**-------------------------------------------------------------------------*/
char *get_workingdirectory(void )

{

    //char buffer[64*FILENAME_MAX];
	char *bf,*pathfile="$WorkingPath";
    //long len;
	long i;
	short a;
    FILE *istream;
    
    istream=fopen(pathfile,"r");    
    if(istream){
		
		i=0;
		a=0;
		do{
			if(i==0){
				bf = (char *) malloc(sizeof(char));
			}else{
				bf = (char *)realloc(bf,(i+1)*sizeof(char));
			}
			bf[i]=fgetc(istream);
			if(bf[i]==10 || bf[i]==-1){
				a=1;
				bf[i]=0;
			}
			i+=1;
		}while(a==0);
		
		fclose(istream);   

    }else{
		
		/*printf("ENTER THE WORKING DIRECTORY PATH:\n");
		scanf("%s",&buffer);
		len=64*FILENAME_MAX;
    	
    	if(len > (64*FILENAME_MAX)){ 
			t_error("Maximum path length exceeded");
		} else {
			bf=(char *)malloc(len*sizeof(char));
		}
		
		strcpy(bf,buffer);	*/
		
		t_error("You have to specify aworking directory when you run the executable");
	}
	return bf;
}


/**-------------------------------------------------------------------------*/

char *join_strings(char *first, char *second)

{

	char *string;
	int len=strlen(first)+strlen(second)+2;
		
	string=(char*)malloc(len*sizeof(char));	
	string=strcpy(string,first);	
	string=strcat(string,second);
	
	return string;
}

/**-------------------------------------------------------------------------*/


char *get_filename(char *working_directory,char *program)

{

    char *name="\0",*fullname="\0",*pathfile="\0",*S=".inpts";
    FILE *istream;
    static STRINGBIN *strg;
    long position=-1;
    short sign=0;
    
	
    pathfile=join_strings(working_directory,join_strings(program,S));
 	
	if(IO_FILES_COUNTER==0){
		istream=fopen(pathfile,"r");    
   		if(istream){
   		   position=simplefind(istream,"1:");
   		   if(position==-1){
   		   		t_error("String not found");
   		  }else{
   		    	strg=read_stringarray(istream, NOPRINT);
   		    	if(strg==NULL){
					name=(char *)malloc(FILENAME_MAX*sizeof(char));
					if(!name) t_error("There was no allocation space");
					sign++;
					scanf("%s",name);

   		    	}else{
   		    		/* print_stringbin_elements(s,80); */
    				IO_FILES_COUNTER++;
    				name=(strg->co[IO_FILES_COUNTER]+1);
    				WORKING_POSITION=ftell(istream);
    				t_fclose(istream);
    			}

	      }

 	    } else {
			name=(char *)malloc(FILENAME_MAX*sizeof(char));
			if(!name) t_error("There was no allocation space");
			sign++;
			scanf("%s",name);
 	
 		}

	} else if(IO_FILES_COUNTER < strg->index->nh){
			IO_FILES_COUNTER++;
			name=(strg->co[IO_FILES_COUNTER]+1);   
			
	} else {
			name=(char *)malloc(FILENAME_MAX*sizeof(char));
			if(!name) t_error("There was no allocation space");
			sign++;
			scanf("%s",name);
	}
	
	fullname=join_strings(working_directory,name);
	printf("%s\n",fullname); 
       if(sign==1) free(name);
       sign=0;
	   
	return fullname;

}


/**-------------------------------------------------------------------------*/

double get_parameter(char *working_directory,char *program)

{

//char *fullname="\0";
char *pathfile="\0",*S=".inpts";
FILE *istream;
static DOUBLEVECTOR *s;
long position=-1;
//short sign=0;
double parameter=0;

pathfile=join_strings(working_directory,join_strings(program,S));

if(IO_PARMS_COUNTER==0){
	istream=fopen(pathfile,"r");

	if(istream){
		position=simplefind(istream,"2:");
	   	if(position==-1){
	   		t_error("String not found");
		}else{
	    	s=read_doublearray(istream,NOPRINT);
			IO_PARMS_COUNTER++;
			parameter=s->co[IO_PARMS_COUNTER];
			WORKING_POSITION=ftell(istream);
			t_fclose(istream);
   		}
    } else {
		scanf(" %lf",&parameter);
	}
} else if(IO_PARMS_COUNTER < s->nh){
	IO_PARMS_COUNTER++;
	parameter=s->co[IO_PARMS_COUNTER];
} else {
	scanf(" %lf",&parameter);
}
printf("%f\n",parameter); 
return parameter;

}


/**-------------------------------------------------------------------------*/

char *get_strings(char *working_directory,char *program)

{

    char *name="\0",*pathfile="\0",*S=".inpts";
	//char *fullname="\0";
    FILE *istream;
    static STRINGBIN *s;
    long position=-1;
    short sign=0;
    
    pathfile=join_strings(working_directory,join_strings(program,S));
    if(IO_FILES_COUNTER==0){
    	        istream=fopen(pathfile,"r");    
   		if(istream){
   		   position=simplefind(istream,"1:");
   		   if(position==-1){
   		   		t_error("String not found");
   		  }else{
   		    	s=read_stringarray(istream, NOPRINT);
   		    	if(s==NULL){
					name=(char *)malloc(FILENAME_MAX*sizeof(char));
					if(!name) t_error("There was no allocation space");
					sign++;
					scanf("%s",name);

   		    	}else{
   		    		/* print_stringbin_elements(s,80); */
    				IO_FILES_COUNTER++;
    				name=(s->co[IO_FILES_COUNTER]+1);
    				WORKING_POSITION=ftell(istream);
    				t_fclose(istream);
    			}

   		}
 	    } else {
			name=(char *)malloc(FILENAME_MAX*sizeof(char));
			if(!name) t_error("There was no allocation space");
			sign++;
			scanf("%s",name);
 	
 		}

	} else if(IO_FILES_COUNTER < s->index->nh){
			IO_FILES_COUNTER++;
			name=(s->co[IO_FILES_COUNTER]+1);   
	} else {
			name=(char *)malloc(FILENAME_MAX*sizeof(char));
			if(!name) t_error("There was no allocation space");
			sign++;
			scanf("%s",name);
	}
	

	printf("%s\n",name); 
	return name;

}

/**-----------------------------------------------------------------------*/

STRINGBIN *read_filenames(char *working_directory,char *program, char *extension, char *position){
	
	STRINGBIN *strg;
	FILE *istream;
	long i,pos=-1;
	
	istream=fopen(join_strings(working_directory,join_strings(program,extension)),"r");
	
	if(istream){
		pos=simplefind(istream,position);
		if(pos==-1) t_error("string not found");
		strg=read_stringarray(istream, NOPRINT);
		fclose(istream);
	}else{
		t_error("File does not exist");
	}
	
	for(i=1;i<=strg->index->nh;i++){
		strg->co[i]=join_strings(working_directory,strg->co[i]+1)-1;
	}
	
	return strg;
}

/**-----------------------------------------------------------------------*/
DOUBLEVECTOR *read_parameters(char *working_directory,char *program, char *extension, char *position){
	
	DOUBLEVECTOR *v;
	FILE *istream;
	long pos=-1;
	
	istream=fopen(join_strings(working_directory,join_strings(program,extension)),"r");
	
	if(istream){
		pos=simplefind(istream,position);
		if(pos==-1) t_error("string not found");
		v=read_doublearray(istream, NOPRINT);
		fclose(istream);
	}else{
		t_error("File does not exist");
	}
	
	return v;
}


/**-----------------------------------------------------------------------*/

char *query_for_label(FILE *inputfile)



{

char ch,ocurl='{';
static char *label=NULL;
long i=0,label_size=0;


		label_size=LABELINCREMENT;

		label=(char *)malloc(label_size*sizeof(char));

		if(label==NULL){

			t_error("Cannot allocate the label");

		}



	 skip_whitespaces(inputfile);

	do{

	  if(i< label_size-1) {
	  		
	  		
	  		label[i]=fgetc(inputfile);
	  		
	  		if(label[i]==EOF){
	  	
	  			t_error("Unespected End Of File");
	  		
	  		}
	
	 		i++; 
	 	

	  } else {
	    
	  		if(label_size==MAXLABELSIZE){
	  		
	  			printf("\nWarning:: the name of the stored data exceeded the LABELSIZE\n");
	  			
	  			printf("Warning:: flushing the rest of the name\n");
	  		
	  			do{
	  			
	  				ch=fgetc(inputfile);
	  			
	  			}while(ch!=ocurl && ch!=EOF);
	  		
	  			if(ch==EOF){
	  			
	  				t_error("Unespected End Of File");
	  			
	  			}
	  			
	  			label[i-1]=ch;
	  			
	  		}else if(label_size + LABELINCREMENT > MAXLABELSIZE){
	    	    		
				label_size=MAXLABELSIZE;
				
	  			label=(char *)realloc(label,label_size*sizeof(char));

	  			if(label==NULL){

	  				t_error("Cannot expand the label");


	  			}
	    	    		
	    	    		
	    	
	    	} else {

				label_size+=LABELINCREMENT;
				
	  			label=(char *)realloc(label,label_size*sizeof(char));

	  			if(label==NULL){

	  				t_error("Cannot expand the label");


	  			}
			}
	  }

	} while(label[i-1]!=EOF && label[i-1]!=ocurl);
	
	ungetc(label[i-1],inputfile);
	label[i]=label[i-1]='\0';
	i-=2;
 
   while(isspace(label[i])){
        label[i]='\0';
    	i--;
    }
    
    i++;    

    
	if(i<=0){
	
	label="NOLABEL";
		 
	}

	return label;




}

/**-----------------------------------------------------------------------*/

char *get_phrase(FILE *inputfile,const char separator)



{



static char *label=NULL,ch;

long i=0,label_size=0;


		label_size=LABELINCREMENT;

		label=(char *)malloc(label_size*sizeof(char));

		if(label==NULL){

			t_error("Cannot allocate the label");

		}



	 skip_whitespaces(inputfile);

	do{

	  if(i< label_size-1) {

	  		label[i]=fgetc(inputfile);

	  		if(label[i]==EOF){
	  	
	  			t_error("Unespected End Of File");
	  		
	  		}
	
	 		i++;
	 	

	  } else {

	  		if(label_size==MAXLABELSIZE){
	  		
	  			printf("\nWarning:: the name of the stored data exceeded the LABELSIZE\n");
	  			
	  			printf("Warning:: flushing the rest of the name\n");
	  		
	  			do{
	  			
	  				ch=fgetc(inputfile);
	  			
	  			}while(ch!=EOF && ch!=separator);
	  		
	  			if(ch==EOF){
	  			
	  				t_error("Unespected End Of File");
	  			
	  			}
	  			
	  			label[i-1]=ch;

	  			
	  		}else if(label_size + LABELINCREMENT > MAXLABELSIZE){
	    	    		
				label_size=MAXLABELSIZE;
				
	  			label=(char *)realloc(label,label_size*sizeof(char));

	  			if(label==NULL){

	  				t_error("Cannot expand the label");


	  			}
	    	    		
	    	    		
	    	
	    	} else {

				label_size+=LABELINCREMENT;
				
	  			label=(char *)realloc(label,label_size*sizeof(char));

	  			if(label==NULL){

	  				t_error("Cannot expand the label");


	  			}
			}
	    

	  }

	} while(/*!isspace(label[i-1]) &&*/label[i-1]!=EOF && label[i-1]!=separator);
	
	/* if( giveback ) ungetc(label[i-1],inputfile); */

	label[i-1]='\0';

    i--;
    
	if(i==0){
	
	label="NOLABEL";
		 
	}
		
	
	/* return label; */
	return join_strings(WORKING_DIRECTORY,label);




}

/**-----------------------------------------------------------------------*/

int query_for_token(FILE *inputfile,const char *token)



{



char keyword[64];

long pos,len,i;

	pos=ftell(inputfile);

    skip_whitespaces(inputfile);

	len=strlen(token);
	
	for(i=0;i<=len-1;i++){

	  keyword[i]=fgetc(inputfile);
	 /* printf("%c",keyword[i]); */

	}

	 keyword[len]='\0';

	if(strcmp(token,keyword)!=0){

		 fseek(inputfile,pos,SEEK_SET);

		return 0;

	}else{

		 return 1;

	}


}


/**-----------------------------------------------------------------------*/





int read_comment(FILE *inputfile,int opt,long maxbuffersize,short print)



{





char *buffer=NULL;



long int buffer_index=0, buffer_size=0;





if(inputfile==NULL){

	

	t_error("You tried to read from a closed file ");



} 



	

buffer_size=BUFFERINCREMENT;



buffer=(char *)malloc(buffer_size*sizeof(char));



if(buffer==NULL){



	t_error("Cannot allocate the buffer");



}





buffer[buffer_size-1]='\0';
buffer[0]=buffer[1]=buffer[2]='\0';


if(opt==1){



	buffer_index=iscomment(buffer,inputfile);



}else{

		

	buffer[0]='/'; buffer[1]='*'; buffer[2]='*';



	buffer_index=3;



}



if (buffer_index==0){
    
    printf("\nWarning:: no comment found\n");

	return 0;
}



do{



	

	if(buffer_index < buffer_size-2){





		buffer[buffer_index]=fgetc(inputfile);



		if(buffer[buffer_index]=='*'){



			buffer_index++;



			buffer[buffer_index]=fgetc(inputfile);



			if(buffer[buffer_index]=='/'){



				buffer_index++;



				buffer[buffer_index]='\0';



				break;



				} else {



					buffer_index++;



				}



		} else {



			++buffer_index;

	

		}



	} else {

	

		if(buffer_size==maxbuffersize){

		



			printf("\nWarning::A very long comment has exceeded the maximum  buffer size\n");



			printf("\nWarnin::Flushing the rest of the comment \n" );



			buffer[buffer_size-8]=' ';

		

			buffer[buffer_size-7]='.';

		

			buffer[buffer_size-6]='.';

		

			buffer[buffer_size-5]='.';

		

			buffer[buffer_size-4]=' ';

		

			do{

				do{

			

					buffer[buffer_size-3]=getc(inputfile);

			

				}while(buffer[buffer_size-3]!='*');

			

				buffer[buffer_size-2]=getc(inputfile);

									

			}while(buffer[buffer_size-2] !='/');

			

			if(buffer[buffer_size-2]==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

		} else {

						

			break;

			

		}

		

		} else if(buffer_size + BUFFERINCREMENT > maxbuffersize){



			buffer_size=maxbuffersize;



			buffer=(char *)realloc(buffer,buffer_size*sizeof(char));



			if(buffer==NULL){



				t_error("Cannot expand the buffer");



			}

			

			buffer[buffer_size-1]='\0';

			

		} else {

		

			buffer_size+=BUFFERINCREMENT;



			buffer=(char *)realloc(buffer,buffer_size*sizeof(char));



			if(buffer==NULL){



				t_error("Cannot expand the buffer");



			}

			

			buffer[buffer_size-1]='\0';



		}

	}



} while( buffer[buffer_index-1]!=EOF );



		

if(buffer[buffer_index-1]==EOF ){



		t_error("Unespected end of file: comment closed by End of file");



}





if(print) printf("\n%s\n",buffer); 



free(buffer);



return buffer_index; 

	



}






/**-----------------------------------------------------------------------*/





char  *readandstore_comment(FILE *inputfile,int opt,long maxbuffersize)


{





char *buffer=NULL;



long int buffer_index=0, buffer_size=0;





if(inputfile==NULL){

	

	t_error("You tried to read from a closed file ");



} 



	

buffer_size=BUFFERINCREMENT;



buffer=(char *)malloc(buffer_size*sizeof(char));



if(buffer==NULL){



	t_error("Cannot allocate the buffer");



}





buffer[buffer_size-1]='\0';



if(opt==1){



	buffer_index=iscomment(buffer,inputfile);



}else{

		

	buffer[0]='/'; buffer[1]='*'; buffer[2]='*';



	buffer_index=3;



}







do{



	

	if(buffer_index < buffer_size-2){





		buffer[buffer_index]=fgetc(inputfile);



		if(buffer[buffer_index]=='*'){



			buffer_index++;



			buffer[buffer_index]=fgetc(inputfile);



			if(buffer[buffer_index]=='/'){



				buffer_index++;



				buffer[buffer_index]='\0';



				break;



				} else {



					buffer_index++;



				}



		} else {



			++buffer_index;

	

		}



	} else {

	

		if(buffer_size==maxbuffersize){

		



			printf("\nWarning::A very long comment has exceeded the maximum  buffer size\n");



			printf("\nWarnin::Flushing the rest of the comment \n" );



			buffer[buffer_size-8]=' ';

		

			buffer[buffer_size-7]='.';

		

			buffer[buffer_size-6]='.';

		

			buffer[buffer_size-5]='.';

		

			buffer[buffer_size-4]=' ';

		

			do{

				do{

			

					buffer[buffer_size-3]=getc(inputfile);

			

				}while(buffer[buffer_size-3]!='*');

			

				buffer[buffer_size-2]=getc(inputfile);

									

			}while(buffer[buffer_size-2] !='/');

			

			if(buffer[buffer_size-2]==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

		} else {

						

			break;

			

		}

		

		} else if(buffer_size + BUFFERINCREMENT > maxbuffersize){



			buffer_size=maxbuffersize;



			buffer=(char *)realloc(buffer,buffer_size*sizeof(char));



			if(buffer==NULL){



				t_error("Cannot expand the buffer");



			}

			

			buffer[buffer_size-1]='\0';

			

		} else {

		

			buffer_size+=BUFFERINCREMENT;



			buffer=(char *)realloc(buffer,buffer_size*sizeof(char));



			if(buffer==NULL){



				t_error("Cannot expand the buffer");



			}

			

			buffer[buffer_size-1]='\0';



		}

	}



} while( buffer[buffer_index-1]!=EOF );



		

if(buffer[buffer_index-1]==EOF ){



		t_error("Unespected end of file: comment closed by End of file");



}


return buffer;

	



}






/**-----------------------------------------------------------------------*/







char *read_buffer_from_stdio(long maxbuffersize)





{



char *bffr=NULL;

long int buffer_index=0, buffer_size=0,h;

if(BUFFERINCREMENT > maxbuffersize){
	printf("\nWarning::BUFFERINCREMENT is larger than MAXBUFFERSIZE\n");
	buffer_size=maxbuffersize;
}else{
	buffer_size=BUFFERINCREMENT;
}


bffr=(char *)malloc(buffer_size*sizeof(char));

if(bffr==NULL){
	t_error("Cannot allocate the buffer");
}



bffr[buffer_size-1]='\0';
bffr[0]='\0';
bffr[1]='\0';
bffr[2]='\0';




h=1;



while(!(bffr[0]=='/' && bffr[1]=='*' && bffr[2]=='*')){

	if(h>1){

		printf("TO START A COMMENT DIGIT '/**' OTHERWISE 'N' \n");		



	}



	do{	

		bffr[0]=getchar();

		if(bffr[0]=='N') {

			bffr[0]='\0';

			return bffr;

		}


	}while(bffr[0]!='/');

		bffr[1]=getchar();

		bffr[2]=getchar();

		h++;


}

buffer_index=0;



do{



	if(buffer_index < buffer_size-1){

		bffr[buffer_index]=getchar();

		buffer_index++;

	} else {



	if(buffer_size==maxbuffersize){

		printf("\nWarning::A very long comment has exceeded the maximum  buffer size\n");

		printf("Warning::Closing the comment\n");

		bffr[buffer_index-2]='\n';

		bffr[buffer_index-3]=' ';

		bffr[buffer_index-4]='.';

		bffr[buffer_index-5]='.';

		bffr[buffer_index-6]='.';

		bffr[buffer_index-7]=' ';

            break;

} else if(buffer_size + BUFFERINCREMENT > maxbuffersize){



	buffer_size=maxbuffersize;

	bffr=(char *)realloc(bffr,buffer_size*sizeof(char));

	if(bffr==NULL){

	t_error("Cannot expand the buffer");



}



	bffr[buffer_size-1]='\0';



} else {



	buffer_size+=BUFFERINCREMENT;

	bffr=(char *)realloc(bffr,buffer_size*sizeof(char));

	if(bffr==NULL){

	t_error("Cannot expand the buffer");



	}

	

	bffr[buffer_size-1]='\0';

}



}



} while( !(bffr[buffer_index-1]=='\n' && bffr[buffer_index-2]=='/' && bffr[buffer_index-3]=='*') );



bffr[buffer_index-3]='\0';



return bffr; 



}





/**-----------------------------------------------------------------------*/

long read_index(FILE *inputfile,short print)

{

char buffer[FILENAME_MAX],curl,ocurl='{',ccurl='}';
long buffer_index=0,pos=0;
//long blocksnumber=0;
long blocks;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
}

while(iscomment(buffer,inputfile)){
	/* read_comment(inputfile,0,MAXBUFFERSIZE,NOPRINT); */
	read_comment(inputfile,0,MAXBUFFERSIZE,print); 
}
pos=ftell(inputfile);
skip_whitespaces(inputfile);

if(query_for_token(inputfile,"index")){

	skip_whitespaces(inputfile);
	curl=fgetc(inputfile);
	
	if(curl!=ocurl) {
		t_error("A non expected character found");
	}

	fscanf(inputfile,"%ld",&blocks);
  	skip_whitespaces(inputfile);

    buffer_index=0;
    buffer[buffer_index]=fgetc(inputfile);
    buffer_index++;
    
    if(buffer[buffer_index-1]==','){
    	buffer_index--;
        do{  
          	buffer[buffer_index]=fgetc(inputfile);
            buffer_index++;
        }while(buffer_index < FILENAME_MAX &&  buffer[buffer_index-1]!=ccurl);
   	}
   	
	if(buffer[buffer_index-1]!=ccurl){
		t_error("Missing closing curl");
	} else {
      	buffer[buffer_index-1]='\0';   
      	if(print) printf("\nTHIS FILE TYPE %s\n\n",buffer);
   	}
}else{
	printf("No blocks specified. One asci block assumed. \n");
    blocks=1;
}

return blocks;

}

/**-----------------------------------------------------------------------*/



int header_scan(FILE* inputfile,HEADER * h)
{
	char buffer[64],semicolon, keyword[MAX_KEYWORD_LENGTH+1],ch;
	char ocurl='{';//ccurl='}';
	//short y=-1,yy=-1;
	long i,j,len,pos;
	extern t_keywords T_KEYWORDS;
	if((h)==NULL){
		t_error("This header was not allocated");
	}
	while(iscomment(buffer,inputfile)){
		read_comment(inputfile,0,MAXBUFFERSIZE,NOPRINT);
	}
	/*
	 read_comment(inputfile,1,MAXBUFFERSIZE,PRINT);
	 */
	pos=ftell(inputfile);
	fscanf(inputfile,"%ld %c",&(h->number),&semicolon);
	if(semicolon!=':'){
		h->number=0;
		fseek(inputfile,pos,SEEK_SET);
	} else {
		pos=ftell(inputfile);
	}
	/* scans the gender, the type and the category of the variable and fills the
	 appropriate subvariable of the header */
	skip_whitespaces(inputfile);
	/* Looking for the gender */
	i=0;
	ch=fgetc(inputfile);
	
	/*printf("^^^%c^^^\n",ch);*/
	
	do{
		keyword[i]=ch;
		ch=fgetc(inputfile);
		i++;
	}while(!isspace(ch) && ch!=ocurl);
	keyword[i]='\0';
	/*  len=strtol(T_KEYWORDS.gender[0],tmp , 10);
	 if(!len){
	 t_error("Error in reading the number of keyword genders");
	 }
	 */
	len=2;
	j=0;
	do{
		j++;
	} while(strcmp(T_KEYWORDS.gender[j],keyword)!=0 && j <=len);
	if(j > len){
		h->gender=1;
	}else{
		h->gender=j;
		if(ch!=ocurl) i=0;
		else i=-1;
	}
	/* Looking for the type */
	if(i==0){
		pos=ftell(inputfile);
		skip_whitespaces(inputfile);
		ch=fgetc(inputfile);
		do{
			keyword[i]=ch;
			ch=fgetc(inputfile);
			i++;
		}while(!isspace(ch) && ch!=ocurl );
		keyword[i]='\0';
	}
	/* len=strtol(T_KEYWORDS.type[0],tmp,10);\
	 if(!len){
	 t_error("Error in reading the number of keywords types");
	 }
	 */
	len=7;
	j=0;
	do{
		j++;
	} while(strcmp(T_KEYWORDS.type[j],keyword)!=0 && j <=len);
	if(j > len){
		h->type=5;/* nel caso la keyword NON sia tra quelle ammesse, viene posta per default a 5 */
	}else{
		h->type=j;
		if(ch!=ocurl) i=0;
		else i=-1;
	}
	/* Looking for the category */
	if(i==0){
		pos=ftell(inputfile);
		skip_whitespaces(inputfile);
		ch=fgetc(inputfile);
		do{
			keyword[i]=ch;
			ch=fgetc(inputfile);
			i++;
		}while(!isspace(ch) && ch!=ocurl);
		keyword[i]='\0';
	}
	/* printf("keyword=%s+\n",keyword); */
	/* len=strtol(T_KEYWORDS.category[0],tmp,10);
	 if(!len){
	 t_error("Error in reading the number of keywords types");
	 }
	 */
	len=5;
	j=0;
	do{
		j++;
		/*printf("%s+%s+\n",T_KEYWORDS.category[j],keyword);*/
	} while(strcmp(T_KEYWORDS.category[j],keyword)!=0 && j <=len);
	if(j > len){
		h->category=3;/* nel caso la keyword NON sia tra quelle ammesse, viene posta per default a 3 */
	}else{
		h->category=j;
		if(ch!=ocurl) i=0;
		else i=-1;
	}
	if(i > 0){
		fseek(inputfile,pos,SEEK_SET);
	} else if( i < 0 ){
		ungetc(ch,inputfile);
	}
	h->name=query_for_label(inputfile);

	return OK;
}
/**-----------------------------------------------------------------------*/

long headercmp(HEADER *a,HEADER *b)



{



short i,count=0;



if(a==NULL || b==NULL){



	t_error("At least one of the headers was not allocated");

	

}



	if(a->gender!=b->gender){

		count++;

	}

	if(a->type!=b->type){

		 count+=10;

	}

	if(a->category!=b->category) { 

		count+=100;

	}

	if(strcmp(a->name,b->name)!=0) {

		count+=1000;

	}		

    if(a->category!=1 && a->category!=4){

		for(i=1;i<=a->dimensions[0];i++){

			if(a->dimensions[i]!=b->dimensions[i]){

				count+=10000;

			}

		}

	}	



	return count;

}



/**-----------------------------------------------------------------------*/





/**-----------------------------------------------------------------------*/



LONGVECTOR *read_longarray(FILE *inputfile,short print)







{







char curl,ch;

const char ocurl='{', ccurl='}';

//short ind=1;

HEADER h;

	//long buffer_index=0;
	long buffer_size=0,i, jj;
	//long blocksnumber=0;
long *numbers;

LONGVECTOR *vec=NULL;





if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} 




/*Scanning the first part of the array */



header_scan(inputfile,&h);



if(print==1){

	print_header(&h);

}

/*Scan for the data in array */

   if((h.category!=1 || h.type!=4) && h.category!=4) 
              printf("\nWarning::the data being read are not stored as an array of long\n");


skip_whitespaces(inputfile);

            

curl=fgetc(inputfile);



if(curl!=ocurl) {



	t_error("A non expected character found");



} else {



	i=0;

	

	skip_whitespaces(inputfile);

	

	buffer_size=BUFFERINCREMENT;

	

	numbers=(long *)malloc(buffer_size*sizeof(long));

    

    if(numbers==NULL){



			t_error("Cannot allocate the buffer");



	}



}



do{

	

	if(i< buffer_size) {

	  



	  		fscanf(inputfile,"%ld ",&numbers[i]);

	  		

	  		ch=fgetc(inputfile);



	  		if(ch==EOF){

	  	

	  			t_error("Unespected End Of File");

	  		

	  		}

	

	 		i++;





	} else {



		if(i==MAXBUFFERSIZE){

			

			printf(" A very long array has exceeded the maximum  buffer size\n");



			printf(" Flushing the rest of the array \n" );

		

			do{

			

					ch=getc(inputfile);

			

			}while(ch!=ccurl && ch!=EOF);

			

			

			if(ch==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

			} else {

						

				break;

			

			}

		

		} else if(buffer_size + BUFFERINCREMENT > MAXBUFFERSIZE){



			buffer_size=MAXBUFFERSIZE;



			numbers=(long *)realloc(numbers,buffer_size*sizeof(long));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;

			

		} else {

		

			buffer_size+=BUFFERINCREMENT;



			numbers=(long *)realloc(numbers,buffer_size*sizeof(long));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;



		}





	}

}while(ch !=EOF && ch!=ccurl );



    

if(i==0){

	

	t_error("Empty array body was encountered");

		 

}





vec=new_longvector(i);



for(jj=0;jj<i;jj++){

    vec->co[jj+1]=numbers[jj]; 

}



free(numbers);





return vec;





}



 /**-----------------------------------------------------------------------*/



FLOATVECTOR *read_floatarray(FILE *inputfile,short print)







{







char curl,ch;

const char ocurl='{', ccurl='}';

//short ind=1;

	//long buffer_index=0;
	long buffer_size=0,i, jj;
//long blocksnumber=0;
float *numbers;

FLOATVECTOR *vec=NULL;

HEADER h;




if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}





/*Scanning the first part of the array */

/* printf("reading a floatarray****\n"); */

header_scan(inputfile,&h);
/*
printf("nnnnn\n");
print_header(&h);
printf("mmmmm\n");
scanf("%c ",&ch);
*/
if(print==1){

	print_header(&h);

}

	if(h.category!=1 || h.type!=5) 
              printf("\nWarning::the data being read are not stored as an array of float\n");
 

/*Scan for the data in array */



skip_whitespaces(inputfile);

            

curl=fgetc(inputfile);



if(curl!=ocurl) {



	t_error("A non expected character found");



} else {



	i=0;

	

	skip_whitespaces(inputfile);

	

	buffer_size=BUFFERINCREMENT;

	

	numbers=(float *)malloc(buffer_size*sizeof(float));

    

    if(numbers==NULL){



			t_error("Cannot allocate the buffer");



	}



}



do{

	

	if(i< buffer_size) {

	  



	  		fscanf(inputfile,"%f ",&numbers[i]);

	  		

	  		ch=fgetc(inputfile);


	  		if(ch==EOF){

	  	

	  			t_error("Unespected End Of File");

	  		

	  		}

	

	 		i++;





	} else {



		if(i==MAXBUFFERSIZE){

			

			printf(" A very long array has exceeded the maximum  buffer size\n");



			printf(" Flushing the rest of the array \n" );

		

			do{

			

					ch=getc(inputfile);

			

			}while(ch!=ccurl && ch!=EOF);

			

			

			if(ch==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

			} else {

						

				break;

			

			}

		

		} else if(buffer_size + BUFFERINCREMENT > MAXBUFFERSIZE){



			buffer_size=MAXBUFFERSIZE;



			numbers=(float *)realloc(numbers,buffer_size*sizeof(float));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;

			

		} else {

		

			buffer_size+=BUFFERINCREMENT;



			numbers=(float *)realloc(numbers,buffer_size*sizeof(float));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;



		}





	}

}while(ch !=EOF && ch!=ccurl );



    

if(i==0){

	

	t_error("Empty array body was encountered");

		 

}





vec=new_floatvector(i);



for(jj=0;jj<i;jj++){

    vec->co[jj+1]=numbers[jj]; 

}



free(numbers);





return vec;





}

 /**-----------------------------------------------------------------------*/



DOUBLEVECTOR *read_doublearray(FILE *inputfile,short print)







{







char curl,ch;

const char ocurl='{', ccurl='}';

//short ind=1;

	//long buffer_index=0;
	long buffer_size=0,i, jj;
//long blocksnumber=0;
double *numbers;

DOUBLEVECTOR *vec=NULL;

HEADER h;



if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} 






/*Scanning the first part of the array */



header_scan(inputfile,&h);


if(print==1){

	print_header(&h);

}


if(h.category!=1 || h.type!=6) 
              printf("\nWarning::the data being read are not stored as an array of double\n");


/*Scan for the data in array */



skip_whitespaces(inputfile);

            

curl=fgetc(inputfile);



if(curl!=ocurl) {



	t_error("A non expected character found");



} else {



	i=0;

	

	skip_whitespaces(inputfile);

	

	buffer_size=BUFFERINCREMENT;

	

	numbers=(double *)malloc(buffer_size*sizeof(double));

    

    if(numbers==NULL){



			t_error("Cannot allocate the buffer");



	}



}



do{

	

	if(i< buffer_size) {

	  



	  		fscanf(inputfile,"%lf ",&numbers[i]);

	  		

	  		ch=fgetc(inputfile);



	  		if(ch==EOF){

	  	

	  			t_error("Unespected End Of File");

	  		

	  		}

	

	 		i++;





	} else {



		if(i==MAXBUFFERSIZE){

			

			printf(" A very long array has exceeded the maximum  buffer size\n");



			printf(" Flushing the rest of the array \n" );

		

			do{

			

					ch=getc(inputfile);

			

			}while(ch!=ccurl && ch!=EOF);

			

			

			if(ch==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

			} else {

						

				break;

			

			}

		

		} else if(buffer_size + BUFFERINCREMENT > MAXBUFFERSIZE){



			buffer_size=MAXBUFFERSIZE;



			numbers=(double *)realloc(numbers,buffer_size*sizeof(double));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;

			

		} else {

		

			buffer_size+=BUFFERINCREMENT;



			numbers=(double *)realloc(numbers,buffer_size*sizeof(double));



			if(numbers==NULL){



				t_error("Cannot expand the buffer");



			}

			

			numbers[buffer_size-1]=0;



		}





	}

}while(ch !=EOF && ch!=ccurl );



    

if(i==0){

	

	t_error("Empty array body was encountered");

		 

}





vec=new_doublevector(i);



for(jj=0;jj<i;jj++){

    vec->co[jj+1]=numbers[jj]; 

}



free(numbers);





return vec;





}



 /**-----------------------------------------------------------------------*/

SHORTVECTOR *read_shortarray(FILE *inputfile,short print)

{


char curl,ch;
const char ocurl='{', ccurl='}';
//short ind=1;
//long buffer_index=0;
long buffer_size=0,i, jj;
//long blocksnumber=;
short *numbers;
SHORTVECTOR *vec=NULL;
HEADER h;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} 

/*Scanning the first part of the array */
header_scan(inputfile,&h);
if(print==1){
	print_header(&h);
}

if(h.category!=1 || h.type!=2) 
              printf("\nWarning::the data being read are not stored as an array of short\n");

/*Scan for the data in array */
skip_whitespaces(inputfile);
curl=fgetc(inputfile);

if(curl!=ocurl) {
	t_error("A non expected character found");
} else {
	i=0;
	skip_whitespaces(inputfile);
	buffer_size=BUFFERINCREMENT;
	numbers=(short *)malloc(buffer_size*sizeof(short));
    if(numbers==NULL){
			t_error("Cannot allocate the buffer");

	}
}

do{
	if(i< buffer_size) {
			fscanf(inputfile,"%hd ",&numbers[i]);
			ch=fgetc(inputfile);
	  		if(ch==EOF){
	  			t_error("Unespected End Of File");
	  		}
	 		i++;
	} else {
		if(i==MAXBUFFERSIZE){
			printf(" A very long array has exceeded the maximum  buffer size\n");
			printf(" Flushing the rest of the array \n" );
			do{
					ch=getc(inputfile);
			}while(ch!=ccurl && ch!=EOF);
			if(ch==EOF){
				t_error("Unespected end of file: comment closed by End of file");
			} else {
				break;
			}
		} else if(buffer_size + BUFFERINCREMENT > MAXBUFFERSIZE){
			buffer_size=MAXBUFFERSIZE;
			numbers=(short *)realloc(numbers,buffer_size*sizeof(short));
			if(numbers==NULL){
				t_error("Cannot expand the buffer");
			}
			numbers[buffer_size-1]=0;
		} else {
			buffer_size+=BUFFERINCREMENT;
			numbers=(short *)realloc(numbers,buffer_size*sizeof(short));
			if(numbers==NULL){
				t_error("Cannot expand the buffer");
			}
			numbers[buffer_size-1]=0;
		}
	}
}while(ch !=EOF && ch!=ccurl );

if(i==0){
	t_error("Empty array body was encountered");

}


vec=new_shortvector(i);

for(jj=0;jj<i;jj++){
    vec->co[jj+1]=numbers[jj]; 
}

free(numbers);

return vec;

}

/**-----------------------------------------------------------------------*/



STRINGBIN *read_stringarray(FILE *inputfile,short print)





{







char curl;

const char ocurl='{', ccurl='}';

//short ind=1;

long strings_size=0,char_counter=0,strings_counter=0;

	long numbers_size=0;
	//long number_index=0,blocksnumber=0,strings_index=0;

char *strings,ch;

long *numbers,iu,ju,sum=0;


HEADER h;

LONGVECTOR *vec=NULL;

STRINGBIN *ST=NULL;





if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} 





/*Scanning the first part of the array */



header_scan(inputfile,&h);


if(print==1){

	print_header(&h);

}


if(h.category!=1 || h.type!=7) 
              printf("\nWarning::the data being read are not stored as an array of strings\n");


/*Scan for the data in array */



skip_whitespaces(inputfile);

            

curl=fgetc(inputfile);



if(curl!=ocurl) {



	t_error("A non expected character found");



} else {



	char_counter=0;

	

	skip_whitespaces(inputfile);

	

	strings_size=BUFFERINCREMENT;

	numbers_size=BUFFERINCREMENT;

	

	strings=(char *)malloc((strings_size+1)*sizeof(char));

    numbers=(long *)malloc(numbers_size*sizeof(long));

    

    if(numbers==NULL || strings==NULL){



			t_error("Cannot allocate the buffers");



	}



}



strings[strings_size-1]='\0';



do{

	

	if(char_counter< strings_size-1 && strings_counter < numbers_size-1) {

	  

	  			

	  		strings[char_counter]=fgetc(inputfile);
/* printf("%d:%c+",char_counter,strings[char_counter]); */
	  		if(strings[char_counter]==EOF){

	  	

	  			t_error("Unespected End Of File");

	  		

	  		} else if (strings[char_counter]==','){

	  		

	  			while(strings[char_counter-1]==' '){

	  				char_counter--;

	  			}

	  		    

	  		    numbers[strings_counter]=char_counter;

	  		   

	  		   	strings_counter++;

	  		   	  		    

	  		    skip_whitespaces(inputfile);

	  			

	  		} else {

	

	 		char_counter++;

	 		

	 		}

	

	  } else if(strings_counter >= numbers_size-1){

		

				numbers_size+=BUFFERINCREMENT;			  	    

				numbers=(long *)realloc(numbers,numbers_size*sizeof(long));



	  			if(numbers==NULL){

					t_error("Cannot expand the buffer");

				}

		

	 } else if(char_counter==MAXBUFFERSIZE-1){

			

			printf("\nWarning::A very long array has exceeded the maximum  buffer size\n");



			printf("Warning::Flushing the rest of the array \n" );

		

			do{

			

					ch=getc(inputfile);

			

			}while(ch!=ccurl && ch!=EOF);

			

			

			if(ch==EOF){

			

				t_error("Unespected end of file: comment closed by End of file");

			

			} else {

			

				strings[char_counter-1]=ch;

						

				break;

			

			}

		

	} else if(strings_size + BUFFERINCREMENT > MAXBUFFERSIZE){



		strings_size=MAXBUFFERSIZE;



		strings=(char *)realloc(strings,(strings_size+1)*sizeof(char));



		if(strings==NULL){



			t_error("Cannot expand the buffer");



		}

			

		strings[MAXBUFFERSIZE-1]='\0';

			



		} else {

		

			strings_size+=BUFFERINCREMENT;



			strings=(char *)realloc(strings,(strings_size+1)*sizeof(char));



			if(strings==NULL){



				t_error("Cannot expand the buffer");



			}

			

			strings[strings_size-1]='\0';

	}



}while(strings[char_counter-1] !=EOF && strings[char_counter-1]!=ccurl );





/* ungetc(strings[char_counter],inputfile); */



char_counter--; 

   

if(char_counter==0){

	

	/* t_error("Empty array body was encountered"); */
       printf("\nWarning :: Empty array body was encountered\n");
	
	return NULL;	 

} else if( strings[char_counter]==ccurl){



  numbers[strings_counter]=char_counter;

    if(char_counter >1){
		while(char_counter >0 && isspace(strings[char_counter-1])){
			char_counter--;
		}
	}
  strings[char_counter]='\0';



} else {



	t_error("Unespected End of File Encountered");



}







vec=new_longvector(strings_counter+1);



vec->co[1]=numbers[0]+1;



for(iu=1;iu<=strings_counter;iu++){

	vec->co[iu+1]=numbers[iu]-numbers[iu-1]+1;

}



free(numbers); 



/* print_longvector_elements(vec,10); */



ST=new_stringbin(vec); 



/* print_longvector_elements(ST->index,10); */



   for(iu=1;iu<=(ST->index)->nh;iu++){

		

		

		for(ju=1;ju< (ST->index)->co[iu];ju++){



			ST->co[iu][ju]=strings[sum];

			

			sum++;

		}	

	

		ST->co[iu][ju]='\0'; 

	

}



free(vec);



return ST;





}






/**-----------------------------------------------------------------------*/



void justread_longarray(FILE *inputfile,LONGVECTOR *vec,short print)

{

char curl,ch;
const char ocurl='{', ccurl='}';
HEADER h;
long i;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");

} 

if(vec==NULL || vec->isdynamic!=1){

	t_error("This vector was not properly allocated");
}

/*Scanning the first part of the array */
header_scan(inputfile,&h);
if(print==1){
	print_header(&h);

}

if(h.category!=1 || !(h.type==4 || h.type==3 || h.type==2)) 
              printf("\nWarning::the data being read cannot be stored as an array of long\n");

/*Scan for the data in array */
skip_whitespaces(inputfile);
	curl=fgetc(inputfile);
	if(curl!=ocurl) {
		t_error("A non expected character found");
	} 
	for(i=1;i<=vec->nh;i++){
	
		ch=fscanf(inputfile," %ld %c",&(vec->co[i]),&curl);

	}

	if(curl!=ccurl){

		t_error("Closing bracket not found");
	} else if(ch==EOF){
	
		t_error("Unespected End of file");	
      }
}



/**-----------------------------------------------------------------------*/

void justread_floatarray(FILE *inputfile,FLOATVECTOR *vec,short print)

{

char curl,ch;
const char ocurl='{', ccurl='}';
HEADER h;
long i;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} 

if(vec==NULL || vec->isdynamic!=1){
	t_error("This vector was not properly allocated");
}

header_scan(inputfile,&h);
if(print==1){
	print_header(&h);
}

if(h.category!=1 || !(h.type==5 || h.type==6)) 
	printf("\nWarning::the data being read cannot be stored as an array of long\n");

skip_whitespaces(inputfile);
curl=fgetc(inputfile);
if(curl!=ocurl) {
	t_error("A non expected character found");
} 
for(i=1;i<=vec->nh;i++){
	ch=fscanf(inputfile," %f %c",&(vec->co[i]),&curl);
}
if(curl!=ccurl){
	t_error("Closing bracket not found");
} else if(ch==EOF){
	t_error("Unespected End of file");	
}
}

/**-----------------------------------------------------------------------*/

void justread_floatmatrix(FILE *inputfile, FLOATMATRIX *C, char *mode,short print)

{
char ch;
HEADER h;
/*FLOATMATRIX * C;*/
long u; 

if(inputfile==NULL){ 
	t_error("You tried to read from a closed file ");
} else { 
	read_matrixheader(inputfile, &h); 
	if(print==1){ 
		print_header(&h); 
		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]); 
	} 
	if(h.category!=3 || h.type!=5) 
		printf("\nWarning::the data being read are not stored as a matrix of float\n"); 
	/***/ 
	/*C=new_floatmatrix(h.dimensions[1],h.dimensions[2]);*/ 
	skip_whitespaces(inputfile); 
	if(query_for_token(inputfile,"->")) { 
		if (OPENYES==1){ 
			t_error("It is not permitted to have more than one external file open"); 
		} else { 
			OPENYES=1; 
			strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';')); 
			if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){ 
				EXTERNAL_FILE_POSITION =SEEK_SET; 
				printf("\nWarning::The external file could open in a wrong position\n"); 
				EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r"); 
			} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){ 
				EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r"); 
				fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET); 
			} else { 
				EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r"); 
				u=read_index(EXTERNAL_FILE,NOPRINT); 
			} 
			read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER); 
			skip_whitespaces(EXTERNAL_FILE); 
			if(headercmp(&h,&EXTERNAL_HEADER)) { 
				t_error("External data file does not match the request "); 
			} else { 
				read_floatmatrix_elements(EXTERNAL_FILE,C,mode); 
				skip_whitespaces(EXTERNAL_FILE); ch=getc(EXTERNAL_FILE); 
				if(ch == EOF){ 
					EXTERNAL_FILE_POSITION=SEEK_SET; 
					t_fclose(EXTERNAL_FILE); 
					free_header(EXTERNAL_HEADER); 
					OPENYES=0; 
				} else { 
					ungetc(ch,EXTERNAL_FILE);
					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE); 
					strcpy(OLD_NAME,EXTERNAL_FILE_NAME); 
					t_fclose(EXTERNAL_FILE); 
					free_header(EXTERNAL_HEADER); 
					OPENYES=0; 
				} 
			} 
		} 
	} else { 
		read_floatmatrix_elements(inputfile,C,mode);
	}
}
free_header(h);
/*return C;*/
}

/**-----------------------------------------------------------------------*/

void justread_chararray(FILE *inputfile,CHARVECTOR *vec,short print)

{

char curl,ch;
const char ocurl='{', ccurl='}';
HEADER h;
long i;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");

} 

if(vec==NULL || vec->isdynamic!=1){

	t_error("This vector was not properly allocated");
}

header_scan(inputfile,&h);
if(print==1){
	print_header(&h);

}


if(h.category!=1 || h.type!=1) 
              printf("\nWarning::the data being read are not stored as an array of char\n");

skip_whitespaces(inputfile);
	curl=fgetc(inputfile);
	if(curl!=ocurl) {
		t_error("A non expected character found");
	} 
	for(i=1;i<=vec->nh;i++){
	
		ch=fscanf(inputfile," %c %c",&(vec->co[i]),&curl);

	}

	if(curl!=ccurl){

		t_error("Closing bracket not found");
	} else if(ch==EOF){
	
		t_error("Unespected End of file");	
      }
}





 /**-----------------------------------------------------------------------*/



void read_vectorheader(FILE *inputfile,HEADER *h)







{







char curl;

const char ocurl='{', ccurl='}';

//short ind=1;

//long buffer_index=0,buffer_size=0,blocksnumber=0;



if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else if(h==NULL){



	t_error("You addresses a null header ");



} else {





	/*Scanning the first part of the array */



	header_scan(inputfile,h);



	/*Scan for the data in array */



	skip_whitespaces(inputfile);

            

	curl=fgetc(inputfile);



	if(curl!=ocurl) {



		t_error("A non expected character found");



	} else {

    	h->dimensions[0]=1; 

		fscanf(inputfile," %ld ",&(h->dimensions[1]));

    	h->dimensions[2]=0;

    	h->dimensions[3]=0;



	}



	curl=fgetc(inputfile);



	if(curl!=ccurl){

		t_error("Closing delimiter not found");

	}



}



}



 /**-----------------------------------------------------------------------*/



void read_matrixheader(FILE *inputfile,HEADER *h)







{







char curl;

const char ocurl='{', ccurl='}';

//short ind=1;

//long buffer_index=0,buffer_size=0,blocksnumber=0;





/*Scanning the first part of the array */



if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else if(h==NULL){



	t_error("You addressed a null header ");



} else {



	header_scan(inputfile,h);


	skip_whitespaces(inputfile);

            

	curl=fgetc(inputfile);



	if(curl!=ocurl) {



		t_error("A non expected character found");



	} else {

    	h->dimensions[0]=2; 

		fscanf(inputfile," %ld ,  %ld ",&(h->dimensions[1]),&(h->dimensions[2]));

    	h->dimensions[3]=0;



	}



	curl=fgetc(inputfile);



	if(curl!=ccurl){

		t_error("Closing delimiter not found");

	}



}



}



/**-----------------------------------------------------------------------*/

CHARVECTOR *read_charvector(FILE *inputfile,char *mode,short print)



{



char ch;
HEADER h;
CHARVECTOR *C;
long u;


	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {

  
   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimension: %ld \n",h.dimensions[1]);
   
   } 
   
   if(h.category!=2 || h.type!=1) printf("\nWarning::the data being read are not stored as a vector of char\n");
   
   C=new_charvector(h.dimensions[1]);


   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;


   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			


   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }


        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  


         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

  
       			read_charvector_elements(EXTERNAL_FILE,C,mode);  
                                            			        	

        		goto_EOF(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	


    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_charvector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}



/**-----------------------------------------------------------------------*/

SHORTVECTOR* read_shortvector(FILE *inputfile,char *mode, short print)



{



char ch;
HEADER h;
SHORTVECTOR *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {


   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);
   
   } 

   if(h.category!=2 || h.type!=2) printf("\nWarning::the data being read are not stored as an vector of short\n");
 

   C=new_shortvector(h.dimensions[1]);

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_shortvector_elements(EXTERNAL_FILE,C,mode);  

       			        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	

    		

    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_shortvector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}



/**-----------------------------------------------------------------------*/

INTVECTOR* read_intvector(FILE *inputfile,char *mode,short print)



{



char ch;
HEADER h;
INTVECTOR *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {

   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld \n",h.dimensions[1]);
   
   } 
   
   if(h.category!=2 || h.type!=3) printf("\nWarning::the data being read are not stored as a vector of int\n");
 

   C=new_intvector(h.dimensions[1]);
	
   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_intvector_elements(EXTERNAL_FILE,C,mode);  

       			        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	

    		

    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_intvector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}



/**-----------------------------------------------------------------------*/

LONGVECTOR *read_longvector(FILE *inputfile,char *mode, short print)



{



char ch;
HEADER h;
LONGVECTOR *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {

  
   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld\n",h.dimensions[1]);
   
   } 


   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld \n",h.dimensions[1]);
   
   } 
	
   if(h.category!=2 || h.type!=4) 
              printf("\nWarning::the data being read are not stored as a vector of long\n");
 

   C=new_longvector(h.dimensions[1]);

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_longvector_elements(EXTERNAL_FILE,C,mode);  

       			        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	

    		

    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_longvector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}





/**-----------------------------------------------------------------------*/

FLOATVECTOR *read_floatvector(FILE *inputfile,char *mode, short print)



{



char ch;
HEADER h;
FLOATVECTOR *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {


   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld \n",h.dimensions[1]);
   
   } 
   
      if(h.category!=2 || h.type!=5) 
              printf("\nWarning::the data being read are not stored as a vector of float\n");

  
   C=new_floatvector(h.dimensions[1]);

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_floatvector_elements(EXTERNAL_FILE,C,mode);  

       			        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	

    		

    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_floatvector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}







/**-----------------------------------------------------------------------*/

DOUBLEVECTOR *read_doublevector(FILE *inputfile,char *mode, short print)



{



char ch;
HEADER h;
DOUBLEVECTOR *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {

   read_vectorheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld\n",h.dimensions[1]);
   
   } 
   
   if(h.category!=2 || h.type!=6) 
              printf("\nWarning::the data being read are not stored as a vector of double\n");

 
  C=new_doublevector(h.dimensions[1]);

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   		

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_vectorheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



         	skip_whitespaces(EXTERNAL_FILE);



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_doublevector_elements(EXTERNAL_FILE,C,mode);  

       			        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	

    		

    			   

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   				

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_doublevector_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}









/**-----------------------------------------------------------------------*/

SHORTMATRIX* read_shortmatrix(FILE *inputfile, char *mode, short print)



{



char ch;
HEADER h;
SHORTMATRIX *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {

 
    read_matrixheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);
   
   }
   
if(h.category!=3 || h.type!=2) 
              printf("\nWarning::the data being read are not stored as a matrix of short\n");
 
  
   C=new_shortmatrix(h.dimensions[1],h.dimensions[2]);

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));





   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  

        	

        	skip_whitespaces(EXTERNAL_FILE);

        	        	        	

        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_shortmatrix_elements(EXTERNAL_FILE,C,mode);  

        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	   

                    

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;    		

  

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 		

   					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   /*     } 

   	

   	}

   */

   }

   

  } else { /* Not an external file */

   

   		read_shortmatrix_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}







/**-----------------------------------------------------------------------*/

INTMATRIX *read_intmatrix(FILE *inputfile, char *mode,short print)



{



char ch;
HEADER h;
INTMATRIX *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {


 
    read_matrixheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);
   
   } 
   
   if(h.category!=3 || h.type!=3) 
              printf("\nWarning::the data being read are not stored as a matrix of int\n");

  
   C=new_intmatrix(h.dimensions[1],h.dimensions[2]);
   

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));





   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	

        	read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  

        	

        	skip_whitespaces(EXTERNAL_FILE);

        	        	        	

        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        			   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_intmatrix_elements(EXTERNAL_FILE,C,mode);  

        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	   

                    

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;    		

  

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 		

   					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

   }

   

  } else { /* Not an external file */

   

   		read_intmatrix_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}









/**-----------------------------------------------------------------------*/

LONGMATRIX  *read_longmatrix(FILE *inputfile, char *mode, short print)



{



char ch;
HEADER h;
LONGMATRIX *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



} else {

 
    read_matrixheader(inputfile,&h);
   
   if(print==1){
   
   		print_header(&h);
   		
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);
   
   } 
  
if(h.category!=3 || h.type!=4) 
              printf("\nWarning::the data being read are not stored as a matrix of long\n");
  
  
   C=new_longmatrix(h.dimensions[1],h.dimensions[2]);
   

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

   

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;

   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   			EXTERNAL_FILE_POSITION =SEEK_SET;

   			    				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){


   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {


        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	
            u=read_index(EXTERNAL_FILE,NOPRINT); /* & */



        }

            

        	read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  



        	skip_whitespaces(EXTERNAL_FILE);        



        	if(headercmp(&h,&EXTERNAL_HEADER)) {

        		print_header(&h);
        		
        		print_header(&EXTERNAL_HEADER);	   

        		t_error("External data file does not match the request ");	   

        			   

       		} else {	

       			

       			read_longmatrix_elements(EXTERNAL_FILE,C,mode);  

        	

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);



        		

    		if(ch == EOF){	   





  	   				EXTERNAL_FILE_POSITION=SEEK_SET;

    		

    				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;



  				} else {



  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

 		

   					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	  

  	   	   			free_header(EXTERNAL_HEADER);

 	   						

  	   				OPENYES=0;

						

  				} 

        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_longmatrix_elements(inputfile,C,mode);

   		

  }

}

free_header(h);
return C;

}







/**-----------------------------------------------------------------------*/

FLOATMATRIX *read_floatmatrix(FILE *inputfile, char *mode,short print)


{

char ch;
HEADER h;
FLOATMATRIX * C;
long u;	

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} else {
	read_matrixheader(inputfile, &h);   
   	if(print==1){
   		print_header(&h); 	
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);		
   	} 
	if(h.category!=3 || h.type!=5) 
        printf("\nWarning::the data being read are not stored as a matrix of float\n");
     /***/
   	C=new_floatmatrix(h.dimensions[1],h.dimensions[2]);
   	skip_whitespaces(inputfile);
 	if(query_for_token(inputfile,"->")) {
   	if (OPENYES==1){
   		t_error("It is not permitted to have more than one external file open");
   	} else {
   		OPENYES=1;
   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));
   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){
 			 EXTERNAL_FILE_POSITION =SEEK_SET;
   			 printf("\nWarning::The external file could  open in a wrong position\n");
       	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){
  			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);
        } else {
        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	u=read_index(EXTERNAL_FILE,NOPRINT);
        }
    	read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);        		                       
    	skip_whitespaces(EXTERNAL_FILE);
    	if(headercmp(&h,&EXTERNAL_HEADER)) {
    		t_error("External data file does not match the request ");	   
  		} else {
   			read_floatmatrix_elements(EXTERNAL_FILE,C,mode);
    		skip_whitespaces(EXTERNAL_FILE);
    		ch=getc(EXTERNAL_FILE);
   		if(ch == EOF){	   
			EXTERNAL_FILE_POSITION=SEEK_SET;
			t_fclose(EXTERNAL_FILE);
			free_header(EXTERNAL_HEADER);
			OPENYES=0;
		} else {
			ungetc(ch,EXTERNAL_FILE);
			EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);
			strcpy(OLD_NAME,EXTERNAL_FILE_NAME);
			t_fclose(EXTERNAL_FILE);
			free_header(EXTERNAL_HEADER);
			OPENYES=0;
		} 
    } 			
}
  } else { 
   		read_floatmatrix_elements(inputfile,C,mode);
  }
}
free_header(h);
return C;
}


/**-----------------------------------------------------------------------*/

DOUBLEMATRIX *read_doublematrix(FILE *inputfile, char *mode,short print)

{

char ch;
HEADER h;
DOUBLEMATRIX *C;
long u;

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} else {
    read_matrixheader(inputfile,&h); 
	if(print==1){
   		print_header(&h);
   		printf("Dimensions: %ld %ld \n",h.dimensions[1],h.dimensions[2]);
	} 
  
	if(h.category!=3 || h.type!=6) 
	  	printf("\nWarning::the data being read are not stored as a matrix of double\n");
	C=new_doublematrix(h.dimensions[1],h.dimensions[2]);
	skip_whitespaces(inputfile);

	if(query_for_token(inputfile,"->")) {

	   	if (OPENYES==1){
	   		t_error("It is not permitted to have more than one external file open");
	   		
	   	} else {
	  		OPENYES=1;
	  		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));
	  		
	  		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){
	   			 EXTERNAL_FILE_POSITION =SEEK_SET;
	   			 printf("\nWarning::The external file could  open in a wrong position\n");
	        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
	        	 
	   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){
	   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
	        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);
	        	
	        } else {
	        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
	        	u=read_index(EXTERNAL_FILE,NOPRINT);
	        }
	        
	        read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  
	        skip_whitespaces(EXTERNAL_FILE);
	        
	        if(headercmp(&h,&EXTERNAL_HEADER)) {
	        	t_error("External data file does not match the request ");	   
			   
	       	} else {	
	       		read_doublematrix_elements(EXTERNAL_FILE,C,mode);  
	        	skip_whitespaces(EXTERNAL_FILE);
	        	ch=getc(EXTERNAL_FILE);
	    	
	    		if(ch == EOF){	   
	   				EXTERNAL_FILE_POSITION=SEEK_SET;    		
					t_fclose(EXTERNAL_FILE);
					free_header(EXTERNAL_HEADER);
	  	   			OPENYES=0;
	  			
	  			} else {

	  				ungetc(ch,EXTERNAL_FILE);	
	  				EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);
	 				strcpy(OLD_NAME,EXTERNAL_FILE_NAME);
	  	   			t_fclose(EXTERNAL_FILE);
	  	   	   		free_header(EXTERNAL_HEADER);
		   			OPENYES=0;
	 			} 
	        } 	
   		}
	} else { /* Not an external file */
   		read_doublematrix_elements(inputfile,C,mode);
	}

}

free_header(h);
return C;

}


/**-----------------------------------------------------------------------*/

INTBIN *read_intbin(FILE *inputfile, char *mode,short print)



{



char ch;
INTBIN* C;
LONGVECTOR *P;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {

    
	P=read_longarray(inputfile,print);
	
	if(print==1){ 
		print_longvector_elements(P,80);

   }


   C=new_intbin(P);
   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,query_for_label(inputfile));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	
        	EXTERNAL_P=read_longarray(EXTERNAL_FILE,NOPRINT);


        	skip_whitespaces(EXTERNAL_FILE);

        	        	

        	if(longvectorcmp(P,EXTERNAL_P)) {

        			   
        		t_error("External data file does not match the request ");	   

        			   

       		} else {


        		free_longvector(EXTERNAL_P);	

        		read_intbin_elements(EXTERNAL_FILE,C,mode);
        		

        		goto_EOF(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);

        		

        		if(ch == EOF){	   

  

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;



  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  	   						

  				} else {

  				

  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

  					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  						

  				}



        	} 	

        			   

   }

 

  } else { 

   

   		read_intbin_elements(inputfile,C,mode);

   		

  }



}


return C;


}



/**-----------------------------------------------------------------------*/

SHORTBIN * read_shortbin(FILE *inputfile, char *mode,short print)



{



char ch;
LONGVECTOR *P;
SHORTBIN *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {


	P=read_longarray(inputfile,print);

	if(print==1){ 
	
		print_longvector_elements(P,80);

   }
   

   C=new_shortbin(P);
   

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,query_for_label(inputfile));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	
        	EXTERNAL_P=read_longarray(EXTERNAL_FILE,NOPRINT);


        	skip_whitespaces(EXTERNAL_FILE);

        	        	

        	if(longvectorcmp(P,EXTERNAL_P)) {

        			   
        		t_error("External data file does not match the request ");	   

        			   

       		} else {


        		free_longvector(EXTERNAL_P);	

        			

        		read_shortbin_elements(EXTERNAL_FILE,C,mode);

        		

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);

        		

        		if(ch == EOF){	   

  

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;



  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  	   						

  				} else {

  				

  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

  					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  						

  				} 



        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_shortbin_elements(inputfile,C,mode);

   		

  }

}


return C;

}





/**-----------------------------------------------------------------------*/

LONGBIN *read_longbin(FILE *inputfile, char *mode, short print)



{



char ch;
LONGVECTOR *P;
LONGBIN *C;
long u;

	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {


	P=read_longarray(inputfile,print);

	if(print==1){ 

		print_longvector_elements(P,80);

   }


   C=new_longbin(P);
   

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,query_for_label(inputfile));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	
        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	
        	EXTERNAL_P=read_longarray(EXTERNAL_FILE,NOPRINT);


        	skip_whitespaces(EXTERNAL_FILE);

        	        	

        	if(longvectorcmp(P,EXTERNAL_P)) {

        			   
        		t_error("External data file does not match the request ");	   

        			   

       		} else {


        		free_longvector(EXTERNAL_P);	

        			

        		read_longbin_elements(EXTERNAL_FILE,C,mode);

        		

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);

        		

        		if(ch == EOF){	   

  

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;



  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  	   						

  				} else {

  				

  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

  					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  						

  				} 



        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_longbin_elements(inputfile,C,mode);

   		

  }

}


return C;

}





/**-----------------------------------------------------------------------*/

DOUBLEBIN *read_doublebin(FILE *inputfile, char *mode,short print)



{



char ch;
LONGVECTOR *P;
DOUBLEBIN *C;
long u;


	

 if(inputfile==NULL){



	t_error("You tried to read from a closed file ");



}  else {

  	P=read_longarray(inputfile,print);

	if(print==1){ 

		print_longvector_elements(P,80);

   }
   


   C=new_doublebin(P);
 

   skip_whitespaces(inputfile);

 

   if(query_for_token(inputfile,"->")) {

        

   	if (OPENYES==1){



   		t_error("It is not permitted to have more than one external file open");



   	} else {

   		

   		OPENYES=1;



   		strcpy(EXTERNAL_FILE_NAME,query_for_label(inputfile));



   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){

   				

   		

   			 EXTERNAL_FILE_POSITION =SEEK_SET;

   				

   			 printf("\nWarning::The external file could  open in a wrong position\n");



        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

   				

   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){

   			

   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");

        	

        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        		

        } else {

        	
        	

        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        	        	
        	u=read_index(EXTERNAL_FILE,NOPRINT);

        }

        	
        	EXTERNAL_P=read_longarray(EXTERNAL_FILE,NOPRINT);


        	skip_whitespaces(EXTERNAL_FILE);

        	        	

        	if(longvectorcmp(P,EXTERNAL_P)) {

        			   
        		t_error("External data file does not match the request ");	   

        			   

       		} else {


        		free_longvector(EXTERNAL_P);	

        		read_doublebin_elements(EXTERNAL_FILE,C,mode);

        		

        		skip_whitespaces(EXTERNAL_FILE);

        		

        		ch=getc(EXTERNAL_FILE);

        		

        		if(ch == EOF){	   

  

  	   				EXTERNAL_FILE_POSITION=SEEK_SET;



  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  	   						

  				} else {

  				

  					ungetc(ch,EXTERNAL_FILE);

  						

  					EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);

  							

  					strcpy(OLD_NAME,EXTERNAL_FILE_NAME);

  							

  	   				t_fclose(EXTERNAL_FILE);

  	   						

  	   				free_header(EXTERNAL_HEADER);

  	   						

  	   				OPENYES=0;

  						

  				} 



        	} 	

        			   

   }

   

  } else { /* Not an external file */

   

   		read_doublebin_elements(inputfile,C,mode);

   		

  }

}


return C;

}








/**-----------------------------------------------------------------------*/

void write_shortarray_elements(FILE *outputfile,SHORTVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%hd,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%hd}\n",V->co[i]);

	

}



putchar('\n');







}





/**-----------------------------------------------------------------------*/

void write_intarray_elements(FILE *outputfile,INTVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%d,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%d}\n",V->co[i]);

	

}



putchar('\n');







}







/**-----------------------------------------------------------------------*/

void write_longarray_elements(FILE *outputfile,LONGVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%ld,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%ld}\n",V->co[i]);

	

}



putchar('\n');







}





/**-----------------------------------------------------------------------*/

void write_floatarray_elements(FILE *outputfile,FLOATVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%f,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%f}\n",V->co[i]);

	

}



putchar('\n');







}



/**-----------------------------------------------------------------------*/

void write_doublearray_elements(FILE *outputfile,DOUBLEVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%lf,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%lf}\n",V->co[i]);

	

}



putchar('\n');







}





/**-----------------------------------------------------------------------*/

void write_chararray_elements(FILE *outputfile, CHARVECTOR *V, long columns)



{



long i;





putchar(' ');



if (V==NULL || V->co==NULL || V->isdynamic !=1){

	t_error("The vector was not allocated properly");

}else if(V->nl > V->nh ){

	t_error("The vector has no proper dimensions");

} else {

    

    fprintf(outputfile,"{");

    

	for(i=V->nl;i<V->nh;i++){

			fprintf(outputfile,"%c,",V->co[i]);

			if(i%columns==0) putchar('\n');

	}

	

		fprintf(outputfile,"%c}\n",V->co[i]);
}
putchar('\n');
}



/**-----------------------------------------------------------------------*/

void write_stringsarray_elements(FILE *outputfile,STRINGBIN *V)
{
long i;
putchar(' ');
if (V==NULL || V->co==NULL || V->isdynamic !=1){
	t_error("The vector was not allocated properly");
} else {
    fprintf(outputfile,"{");
	for(i=(V->index)->nl;i<(V->index)->nh;i++){
			fprintf(outputfile,"%s,\n",V->co[i]);
	}
		fprintf(outputfile,"%s}\n",V->co[i]);
}
putchar('\n');
}






/**-----------------------------------------------------------------------*/



void print_header(HEADER *H)
{
extern t_keywords T_KEYWORDS;
if(H != NULL){
	        printf("\nBlock number: %ld.\nType:  %s %s %s\nName: %s\n",H->number,  T_KEYWORDS.gender[H->gender], T_KEYWORDS.type[H->type], T_KEYWORDS.category[H->category], H->name);
}else {
	printf("\nWarning::This header was not allocated\n");
}
}



/**-----------------------------------------------------------------------*/

void write_header_header(FILE *outputfile,HEADER *H)
{
long i;
if(outputfile==NULL || H==NULL || H->name==NULL  ){
	t_error("An attempt was made to write on a not opened file or a NULL vector");
}else {
	fprintf(outputfile,"%ld: %s %s %s %s ",H->number,T_KEYWORDS.gender[H->gender],T_KEYWORDS.type[H->type], T_KEYWORDS.category[H->category], H->name);
	 if(H->dimensions[0]!=0 && H->dimensions[0] < 3){
		fprintf(outputfile,"%s",T_KEYWORDS.delimiter[1]);
		for(i=1;i<H->dimensions[0];i++){
				fprintf(outputfile,"%ld,",H->dimensions[i]);	
		}	 
		fprintf(outputfile,"%ld%s\n",H->dimensions[i],T_KEYWORDS.delimiter[2]);		   
	 }
 }
}









/**-----------------------------------------------------------------------*/

void write_turtle(FILE * outputfile,char *creator, char *inputs)
{
if(outputfile==NULL){
	t_error("An attempt was made towrite on a not opened file");
} else if(inputs ==NULL && creator==NULL) {
	fprintf(outputfile,"%3s %21s created on %s at %s %2s \n",\
	 "/**","This_is_a_turtle_file",__DATE__,__TIME__,"*/");
} else if(inputs ==NULL) {
	fprintf(outputfile,"%3s %21s created on %s at %s by %s %2s\n",\
	 "/**","This_is_a_turtle_file",__DATE__,__TIME__,creator,"*/");
} else if(creator==NULL) {
	fprintf(outputfile,"%3s %21s created on %s at %s\n inputs processed: %s %2s\n","/**","This_is_a_turtle_file",__DATE__,__TIME__,inputs,"*/");
} else {
	 fprintf(outputfile,"%3s %21s created on %s at %s  by %s \n inputs processed :%s %s\n",\
	 "/**","This_is_a_turtle_file",__DATE__,__TIME__,creator,inputs,"*/");
}
}



/**-----------------------------------------------------------------------*/

void write_comment(FILE * outputfile,const char *comment, long columns)
{
const char  opening[5]="/** ",closing[4]=" */";
long i,j,len;
if(outputfile==NULL){
	t_error("An attempt was made to write on a not opened file");
} else if(comment!=NULL) {
   len=strlen(comment);
    /* buffer=(char *)malloc((columns+1+NR_END)*sizeof(char)); */
    putc('\n',outputfile);
    fprintf(outputfile,"%s",opening);
	j=5;
	for(i=0;i<len;i++){
		if( j< columns){
		            if(  comment[i]!='\n' && comment[i]!='\t'){
						putc(comment[i],outputfile);
						} else {
						putc(' ',outputfile);
						}
					j++;
			}
		if(j == columns && !isspace(comment[i-1])){
			putc('-',outputfile);
			putc('\n',outputfile); 
			j=0;
		}
	}
	if( j+2 < columns){
		fprintf(outputfile,closing);
    } else {
    	putc('\n',outputfile);
    	fprintf(outputfile,closing);
    }
    putc('\n',outputfile);
} 
}


/*-----------------------------------------------------------------------*/

long read_doubletensor_elements(FILE *input,DOUBLETENSOR *m,char *mode)
/** Read a matrix of double numbers */

{

long tmp=0,count=0;

const char ascii[2]="a",binary[2]="b";
long i,j,k;
if(input==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || (m->isdynamic)!=1){
	t_error("The tensor was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl > m->nch || m->ndl > m->ndh){
	t_error("The tensor has no proper dimensions");
} else if(strcmp(mode,ascii)==0){

for(k=m->ndl;k<=m->ndh;k++){	
	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
				tmp=fscanf(input,"%lf",&(m->co[k][i][j]));
				if(tmp!=EOF){
					count+=tmp;
				}else {
					printf("Error in stored data::Unexpected End of File encountered\n");
					printf(" after position %ld\n",count);
					return -count;
				}
			}
		}
	}
	
}else if(strcmp(mode,binary)==0){
	t_error("Error in reading mode::Mode not yet implemented");  	
} else {
	    t_error("Error in reading mode::Mode not supported");
}

if(count!=(m->nrh-m->nrl+1)*(m->nch-m->ncl+1)*(m->ndh-m->ndl+1)){
	printf("Error in stored data::Stored data number does not match the request");
	return -count;
}else{
	return count;
}

}


/**-----------------------------------------------------------------------*/
long write_doubletensor_elements(FILE * output,DOUBLETENSOR *m,long maxcols)

/* Write a matrix of double numbers to a file */

{

long tmp=0,i,j,k;



if(output==NULL){
	t_error("The input file was not opened properly");
}else if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The tensor was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch || m->ndl > m->ndh){
	t_error("The  tensor has no proper dimensions");
} else {

for(k=m->ndl;k<=m->ndh;k++){	
	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
				tmp=fprintf(output,"%f ",m->co[k][i][j]);
				if(tmp==EOF){
 					printf("Error in storing data::Unespected End of file encountered\n");
					return EOF;
 				}
			} 		
			putc('\n',output);
		    if(j%maxcols==0 && j!= m->ncl) putc('\n',output);
		}
		putc('\n',output);
	}
	putchar('\n');
}

return OK;

}



/**-----------------------------------------------------------------------*/
void print_doubletensor_elements(DOUBLETENSOR *m,long maxcols)
/* Write a matrix of  double to the standard output */

{

long i,j,k;
/* char ch; */

putchar('\n');

if (m==NULL || m->co==NULL || m->isdynamic !=1){
	t_error("The matrix was not allocated properly");
}else if(m->nrl > m->nrh || m->ncl >m->nch ){
	t_error("The matrix has no proper dimensions");
} else {

for(k=m->ndl;k<=m->ndh;k++){			
	for(i=m->nrl;i<=m->nrh;i++){
		for(j=m->ncl;j<=m->nch;j++){
				printf("%f ",m->co[k][i][j]);
		    		if(j%maxcols==0 && j!= m->nch) putchar('\n');
		    	/*	scanf("%c",&ch); */
			} putchar('\n');
		}
		putchar('\n');
	}
	putchar('\n');	
}


}


/**-----------------------------------------------------------------------*/

DOUBLETENSOR *read_doubletensor(FILE *inputfile, char *mode,short print)



{



char ch;
HEADER h;
DOUBLETENSOR *C;
long u;

	

 if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} else {
    read_tensorheader(inputfile,&h);
    if(print==1){
     		print_header(&h);
   		printf("Dimensions: %ld %ld %ld\n",h.dimensions[1],h.dimensions[2],h.dimensions[3]);
   
   } 

  
if(h.category!=5 || h.type!=6) 
    printf("\nWarning::the data being read are not stored as a tensor of double\n");
    C=new_doubletensor(h.dimensions[1],h.dimensions[2],h.dimensions[3]);
   skip_whitespaces(inputfile);
   if(query_for_token(inputfile,"->")) {
   	if (OPENYES==1){
   		t_error("It is not permitted to have more than one external file open");
   	} else {

   		OPENYES=1;
   		strcpy(EXTERNAL_FILE_NAME,get_phrase(inputfile,';'));
   		if(strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION !=SEEK_SET){
   			 EXTERNAL_FILE_POSITION =SEEK_SET;
   			 printf("\nWarning::The external file could  open in a wrong position\n");
	        	 EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
   		} else if(!strcmp(EXTERNAL_FILE_NAME,OLD_NAME) && EXTERNAL_FILE_POSITION!=SEEK_SET){
   			EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
	        	fseek(EXTERNAL_FILE,EXTERNAL_FILE_POSITION,SEEK_SET);

        	} else {
	        	EXTERNAL_FILE=t_fopen(EXTERNAL_FILE_NAME,"r");
        		u=read_index(EXTERNAL_FILE,NOPRINT);
        	}

        	read_matrixheader(EXTERNAL_FILE,&EXTERNAL_HEADER);  
        	skip_whitespaces(EXTERNAL_FILE);
        	if(headercmp(&h,&EXTERNAL_HEADER)) {
        		t_error("External data file does not match the request ");	   
       		} else {	
       			read_doubletensor_elements(EXTERNAL_FILE,C,mode);  
        		skip_whitespaces(EXTERNAL_FILE);
        		ch=getc(EXTERNAL_FILE);
	    		if(ch == EOF){	   
   				EXTERNAL_FILE_POSITION=SEEK_SET;    		
    				t_fclose(EXTERNAL_FILE);
				free_header(EXTERNAL_HEADER);
  	   			OPENYES=0;
			} else {
  				ungetc(ch,EXTERNAL_FILE);
  				EXTERNAL_FILE_POSITION=ftell(EXTERNAL_FILE);
  				strcpy(OLD_NAME,EXTERNAL_FILE_NAME);
  	   			t_fclose(EXTERNAL_FILE);
  	   	   		free_header(EXTERNAL_HEADER);
 	   			OPENYES=0;
				} 
        	} 	
   }

   

  } else { /* Not an external file */
   		read_doubletensor_elements(inputfile,C,mode);
  }

}

free_header(h);
return C;

}



 /**-----------------------------------------------------------------------*/



void read_tensorheader(FILE *inputfile,HEADER *h)
{
char curl;

const char ocurl='{', ccurl='}';

//short ind=1;

//long buffer_index=0,buffer_size=0,blocksnumber=0;

/*Scanning the first part of the array */

if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
} else if(h==NULL){
	t_error("You addressed a null header ");
} else {
	header_scan(inputfile,h);
	skip_whitespaces(inputfile);
	curl=fgetc(inputfile);
	if(curl!=ocurl) {
		t_error("A non expected character found");
	} else {
		h->dimensions[0]=3; 
		fscanf(inputfile," %ld ,  %ld , %ld",&(h->dimensions[1]),&(h->dimensions[2]),&(h->dimensions[3]));
	}

	curl=fgetc(inputfile);
	if(curl!=ccurl){
	t_error("Closing delimiter not found");

	}
}
}









