
/* STATEMENT:

ASCII-GEOtop LIBRARIES

Copyright, 2008 Stefano Endrizzi, Riccardo Rigon

 LICENSE:

 This file is part of ASCII-GEOtop LIBRARIES

 ASCII-GEOtop LIBRARIES is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */


//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <math.h>
//#include <time.h>
//#include <ctype.h>
#include "turtle.h" /* line added by Emanuele Cordano on 1 September 2009 */
#include "tabs.h"

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
char **readline_textarray(FILE *f, long offset){

	long i,icont,j,k;
	char **ch=NULL;

	//printf("offset:%ld\n",offset);

	if(f){
		icont=0;

		do{

			if(icont==0){
				i=0;
				ch = (char **)malloc(sizeof(char*));
			}else if(icont>offset){
				i=icont-offset;
				ch = (char **)realloc(ch,(i+1)*sizeof(char*));
			}else{
				i=0;
				free(ch[i]);
			}

			ch[i] = (char *) malloc(sizeof(char));

			do{
				ch[i][0]=fgetc(f);

				//printf("0: %ld i:%ld icont:%ld\n",ch[i][0],i,icont);
				//stop_execution();

				if(ch[i][0]==47){// ascii character 47='/' which for us means the beginning of a comment
					do{
						ch[i][0]=fgetc(f);
						printf("comment %d\n",ch[i][0]);
						if(ch[i][0]==-1) t_error("ERROR: COMMENT NOT CLOSED");
					}while(ch[i][0]!=47);
					ch[i][0]=fgetc(f);
				}

				//printf("01: %ld i:%ld icont:%ld\n",ch[i][0],i,icont);
				//stop_execution();

				if(i==0 && ch[i][0]==10) ch[i][0]=32;// ascii char=10 means '\n' and 32 means 'Space'
				if(i==0 && ch[i][0]==-1) t_error("ERROR: END OF FILE ENCOUNTERED WHILE READING");
				if(ch[i][0]<=43 || (ch[i][0]>=45 && ch[i][0]<=47) || (ch[i][0]>=58 && ch[i][0]<=64) || (ch[i][0]>=91 && ch[i][0]<=96) || ch[i][0]>=123){
					if(ch[i][0]!=32) t_error("ERROR: NOT ADMITTED COMMENTS NOT PROCEDED BY THE CHARACTER '/' OR COMMENTS NOT CLOSED");
				}
			}while(ch[i][0]<=43 || (ch[i][0]>=45 && ch[i][0]<=47) || (ch[i][0]>=58 && ch[i][0]<=64) || (ch[i][0]>=91 && ch[i][0]<=96) || ch[i][0]>=123);

			j=0;
			//printf("1: i:%ld j:%ld %ld\n",i,j,ch[i][0]);
			//stop_execution();
			if(ch[i][j]!=44 && ch[i][j]!=10 && ch[i][j]!=13 && ch[i][j]!=-1){
				do{
					if( (ch[i][j]>=48 && ch[i][j]<=57) || (ch[i][j]>=65 && ch[i][j]<=90) || (ch[i][j]>=97 && ch[i][j]<=122) ){
						j++;
						//printf("j:%ld\n",j);
						//stop_execution();
						ch[i] = (char *)realloc(ch[i],(j+1)*sizeof(char));
					}
					ch[i][j]=fgetc(f);
					//printf("1: i:%ld j:%ld %ld\n",i,j,ch[i][j]);
					//stop_execution();
				}while(ch[i][j]!=44 && ch[i][j]!=10 && ch[i][j]!=13 && ch[i][j]!=-1);
			}
			if(ch[i][j]==-1) t_error("ERROR: END OF FILE WHILE READING HEADER");

			k=ch[i][j];
			ch[i][j]=0;

			//printf("%s\n",ch[i]);
			//stop_execution();

			icont++;

		}while(k!=10 && k!=13 && k!=-1);

		ch = (char **)realloc(ch,(i+2)*sizeof(char*));
		ch[i+1] = (char *) malloc(sizeof(char));
		ch[i+1][0]=-1;

	}else{
		t_error("ERROR: ATTEMPTING TO READ A CLOSED FILE");
	}

	return(ch);
}



/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
void readline_array(FILE *f, double *a, long offset, long ncol, double ndef, short *endoffile){
	long i,j,k;//,cnt;
	char *ch=NULL;

	if(f){
		*endoffile=0;
		i=0;
		ch = (char *) malloc(sizeof(char));// Matteo 28.9/09
		do{
			//printf("%ld %ld",offset,i);
			//stop_execution();

			//ch = (char *) malloc(sizeof(char));// Matteo 28.9/09
			do{
				ch[0]=fgetc(f);
				//printf("00. %ld i:%ld\n",ch[0],i);
				//stop_execution();
				if(ch[0]==47){	// ascii character 47='/' which for us means a comment
					do{
					//	free(ch); /* commeted by Matteo Dall'Amico on 23/09/09 */
					//	ch = (char *) malloc(sizeof(char));
						ch[0]=fgetc(f);
						//printf("0. %ld i:%ld\n",ch[0],i);
						if(ch[0]==-1) t_error("ERROR: COMMENT NOT CLOSED");
						//stop_execution();
					}while(ch[0]!=47);
				}
				//printf("1. %ld i:%ld\n",ch[0],i);
				// ascii character 10='/n' and 32='Space'
				if(i==0 && ch[0]==10) ch[0]=32;// it means that the first line of the file is a "return"
				if(i==0 && ch[0]==-1) printf("\nWarning: END OF FILE ENCOUNTERED WHILE READING\n");
				if(ch[0]<44 || ch[0]>57){
					if(ch[0]!=32 && ch[0]!=-1) t_error("ERROR: NOT ADMITTED COMMENTS NOT PROCEDED BY THE CHARACTER '/' OR COMMENTS NOT CLOSED");
				}
			}while(ch[0]<-1 || (ch[0]>-1 && ch[0]<10) || (ch[0]>10 && ch[0]<13) || (ch[0]>13 && ch[0]<42) || ch[0]==47 || ch[0]>57);

			j=0;
			if(ch[j]!=44 && ch[j]!=10 && ch[j]!=13 && ch[j]!=-1){
				do{
					j+=1;
					ch = (char *)realloc(ch,(j+1)*sizeof(char));
			//		ch_temp = (char *)realloc(ch,(j+1)*sizeof(char));
			//		free(ch);
			//		ch_temp=ch;

				//	ch_temp=(char *)malloc((j+1)*sizeof(char)); /* modified by Emanuele Cordano on 25 September 2008 */
			//		for(cnt=0;cnt<=j-1;j++){
				//		ch_temp[cnt]=ch[cnt];
				//	}
				//	free(ch);
				//	ch=(char *)malloc((j+1)*sizeof(char));
				//	for(cnt=0;cnt<=j-1;j++){
				//		ch[cnt]=ch_temp[cnt];
				//	}
				//	free(ch_temp);

					//


					//ch_temp=(char *)malloc((j+1)*sizeof(char));

					//free(ch);
					//ch=(char *)malloc((j+1)*sizeof(char));
					//strcpy(ch_temp,ch);
					//free(ch_temp);

					ch[j]=fgetc(f);
					//printf("2. %ld j:%ld\n",ch[j],j);
					//stop_execution();
				}while(ch[j]!=44 && ch[j]!=10 && ch[j]!=13 && ch[j]!=-1);
			}
			if(ch[j]==-1) *endoffile=1;

			if(i>=offset){
				a[i-offset]=decod(ch, j, ndef);
				//printf("J:%ld a(%ld)=%f",j,i-offset,a[i-offset]);
				//stop_execution();
			}
			i++;

		}while(i<ncol+offset && (ch[j]!=10 && ch[j]!=13 && ch[j]!=-1));

		if(ch[j]==44){
		//	free(ch);
		//	ch = (char *) malloc(sizeof(char));
			do{
				ch[0]=fgetc(f);
				//printf("...%ld",ch[0]);
				//stop_execution();
			}while(ch[0]!=10 && ch[j]!=13 && ch[0]!=-1);
			if(ch[j]==-1) *endoffile=1;
		}else{
			for(k=i;k<=ncol+offset-1;k++){
				a[k-offset]=ndef;
				//printf("==a(%ld)=%f",k-offset,a[k-offset]);
				//stop_execution();
			}
		}
	//	free(ch);
	}else{
		t_error("ERROR: ATTEMPTING TO READ A CLOSED FILE");
	}
	free(ch);
}



/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
double decod(char *ch, long n, double ndef){
	short sgn=0;
	long i,ibeg,idec,j;
	double a=0;
	if(ch[0]==42 || ch[0]==44 || ch[0]==-1 || ch[0]==10 || ch[0]==13){
		a=ndef;
	}else{
		i=-1;
		do{
			i+=1;
		}while((ch[i]<46 || ch[i]==47 || ch[i]>57) && (i<n-1));
		if((i==n-1) && (ch[i]<46 || ch[i]==47 || ch[i]>57)){
			a=ndef;
		}else{
			ibeg=i;
			if(ibeg>0){
				if(ch[ibeg-1]==45) sgn=1;
			}
			i-=1;
			do{
				i+=1;
			}while(ch[i]!=46 && (i<n-1));
			if(ch[i]==46){
				idec=i;
			}else{
				idec=n;
			}
			/*if(i!=n-1){
				do{
					i+=1;
					if(ch[i]<48 || ch[i]>57) t_error("ERROR: NUMBER WRITTEN IN A NON-STANDARD FORM");
				}while(i<n-1);
			}*/
			j=0;
			for(i=idec-1;i>=ibeg;i--){
				if(ch[i]>=48 && ch[i]<=57){
					a+=(ch[i]-48)*pow(10.0,(double)(idec-i-1-j));
				}else{
					j+=1;
					if(ch[i]==69) t_error("ERROR:IT IS NOT POSSIBLE TO WRITE NUMBERS IN THE EXPONENTIAL FORMAT");
				}
			}
			j=0;
			for(i=idec+1;i<n;i++){
				if(ch[i]>=48 && ch[i]<=57){
					a+=(ch[i]-48)*pow(10.0,(double)(idec-i+j));
				}else{
					j+=1;
					if(ch[i]==69) t_error("ERROR:IT IS NOT POSSIBLE TO WRITE NUMBERS IN THE EXPONENTIAL FORMAT");
				}
			}

			if(sgn==1) a*=(-1);
		}
	}
	return(a);
}

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
long dim1(double *a){
	long i=0;
	while(a[i]!=end_vector) ++i;
	return(i);
}

long dim1l(long *a){
	long i=0;
	while(a[i]!=end_vector_long) ++i;
	return(i);
}

long dim2(double **a){
	long i=0;
	while(a[i][0]!=end_vector) ++i;
	return(i);
}

long dim_string(char *a){
	long i=0;
	while(a[i]>0) ++i;
	return(i);
}

long dim_vect_strings(char **a){
	long i=0;
	while(a[i][0]>0) ++i;
	return(i);
}

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
double **alloc2(long n, long m){
	/* function that allocates a matrix n*(m+1)
	 * the last column is all populated with "end_vector" numbers
	 * the row n+1 is just a pointer to a double initialized with "end_vector" number
	 * comment: Matteo Dall'Amico, 28/9/09 in Zurich */
	double **a;
	long i;
	a=(double **)malloc( (n+1)*sizeof(double*) );
	for(i=0;i<=n-1;i++){
		a[i]=(double *)malloc( (m+1)*sizeof(double) );
		a[i][m]=end_vector;
	}
	a[n]=(double *)malloc( sizeof(double) );
	a[n][0]=end_vector;
	return(a);
}

double *alloc1(long n){
	double *a;
	long i;
	a=(double *)malloc( (n+1)*sizeof(double) );
	for(i=0;i<=n-1;i++){
		a[i]=0.0;
	}
	a[n]=end_vector;
	return(a);
}

long **alloc_long2(long n){
	long **a,i;
	a=(long **)malloc( (n+1)*sizeof(long*) );
	for(i=0;i<=n;i++){
		a[i]=(long *)malloc( sizeof(long) );
		if(i<n) a[i][0]=0;
		if(i==n) a[i][0]=end_vector_long;
	}
	return(a);
}

long *alloc_long1(long n){
	long *a,i;
	a=(long *)malloc((n+1)*sizeof(long) );
	for(i=0;i<=n;i++){
		if(i<n) a[i]=0;
		if(i==n) a[i]=end_vector_long;
	}
	return(a);
}


/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
short compare_strings(char *a, char *b){
	long i=1;
	short res=0;
	do{
		if(i<=dim_string(a)){
			if(a[i-1]==b[i-1]){
				res=1;
			}else{
				res=0;
			}
			i++;
		}
	}while(res==1 && i<=dim_string(a));
	return res;
}

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
int count_meteo_lines(FILE *f) {
	/* Author: Matteo Dall'Amico, 28.9.09 in Zurich
	* reads a file and counts the lines, without considering the "return" lines
	* Input:
	* f: pointer to the file
	* Output:
	* res: number of lines */
	int c=0;
	int i=0;
	int res=0;
	char ch='\0';// current character
	char ch_prev='\0';// previous character
	if(f){
		while(ch!=EOF) {
			ch_prev=ch;
			ch=fgetc(f);
			if(ch=='\n')  c++;
			if(ch=='\n' && ch_prev=='\n')  i++;
			}
	}else{
			t_error("ERROR: ATTEMPTING TO READ A CLOSED FILE");
		}
	if(ch_prev!='\n') res=c-i+1;
	else res=c-i;
	//printf("ch=%c, ch_prev=%c, c=%d, i=%d, numlines=%d",ch,ch_prev,c,i,res); stop_execution();
	return res; // because 1 line is for the first line of the header
    }


void free_alloc2(double ***_matr)
	/* Author: Matteo Dall'Amico, 28.9.09 in Zurich with the help Dave May
	* deallocates the matrix n*(m+1) allocated by alloc2
	* Input:
	* matr: double pointer to double (the dynamic matrix to deallocate)
	*/
{
	double **matr = *_matr;
	long n=0;// rows
	int i;// counter
	/*for(i=0;i<n;i++){
		for(j=0;j<dim1(matr[0]);j++){
			printf("matr[%d][%d]=%f",i,j,matr[i][j]);
		}
	}
	stop_execution();*/
	while(matr[n][0]!=end_vector) n++;
	//printf("nrows=%ld",n);stop_execution();
	for(i=0;i<n;i++) {
		free(matr[i]);// free each dynamic vector of the matrix
	}
	free(matr);

	*_matr = NULL;
}
