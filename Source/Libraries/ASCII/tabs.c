
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.0 Public - Version "Montebello" - Update 2 (29 April 2010)
 
 Copyright (c), 2010 - Stefano Endrizzi and Riccardo Rigon
 
 This file is part of GEOtop 1.0 Public
 
 GEOtop 1.0 Public is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.0 Public is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community at the following E-mail address: geotopusers@googlegroups.com to which you can subscribe at  http://groups.google.com/group/geotopusers/
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors
 
 */
	
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "tabs.h"

void t_error(char *error_text);

/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------*/
char **readline_textarray(FILE *f, long offset){

	long i,icont,j,k;
	char **ch;
	
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

				if(ch[i][0]==47){	//comment
					do{
						ch[i][0]=fgetc(f);
						printf("comment %d\n",ch[i][0]);
						if(ch[i][0]==-1) t_error("ERROR: COMMENT NOT CLOSED");
					}while(ch[i][0]!=47);
					ch[i][0]=fgetc(f);
				}	

				//printf("01: %ld i:%ld icont:%ld\n",ch[i][0],i,icont);
				//stop_execution();
						
				if(i==0 && ch[i][0]==10) ch[i][0]=32;
				if(i==0 && ch[i][0]==-1) t_error("ERROR: END OF FILE ENCOUNTERED WHILE READING");		
				if(ch[i][0]<=43 || (ch[i][0]>=45 && ch[i][0]<=47) || (ch[i][0]>=58 && ch[i][0]<=64) || (ch[i][0]>=91 && ch[i][0]<=96) || ch[i][0]>=123){
					if(ch[i][0]!=32 && ch[i][0]!=9) t_error("ERROR: NOT ADMITTED COMMENTS NOT PROCEDED BY THE CHARACTER '/' OR COMMENTS NOT CLOSED");
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
	long i,j,k;
	char *ch;
	
	
	if(f){
		*endoffile=0;
		i=0;
		
		do{
			//printf("%ld %ld",offset,i);
			//stop_execution();

			ch = (char *) malloc(sizeof(char));								
			do{
				ch[0]=fgetc(f);
				//printf("00. %ld i:%ld\n",ch[0],i);
				//stop_execution();
				if(ch[0]==47){	//comment
					do{
						ch[0]=fgetc(f);
						//printf("0. %ld i:%ld\n",ch[0],i);
						if(ch[0]==-1) t_error("ERROR: COMMENT NOT CLOSED");
						//stop_execution();
					}while(ch[0]!=47);
				}	
				//printf("1. %ld i:%ld\n",ch[0],i);
						
				if(i==0 && ch[0]==10) ch[0]=32;
				//if(i==0 && ch[0]==-1) printf("\nWarning: END OF FILE ENCOUNTERED WHILE READING\n");	
				if(ch[0]<44 || ch[0]>57){					
					if(ch[0]!=32 && ch[0]!=-1 && ch[0]!=9){
						printf("it happens here\n");
						t_error("ERROR: NOT ADMITTED COMMENTS NOT PROCEDED BY THE CHARACTER '/' OR COMMENTS NOT CLOSED");	
					}
				}
			}while(ch[0]<-1 || (ch[0]>-1 && ch[0]<10) || (ch[0]>10 && ch[0]<13) || (ch[0]>13 && ch[0]<42) || ch[0]==47 || ch[0]>57);
			
			j=0;
			if(ch[j]!=44 && ch[j]!=10 && ch[j]!=13 && ch[j]!=-1){
				do{
					j+=1;
					ch = (char *)realloc(ch,(j+1)*sizeof(char));
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
			free(ch);
			ch = (char *) malloc(sizeof(char));	
			do{
				ch[0]=fgetc(f);
				//printf("...%d",ch[0]);
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
		free(ch);
	}else{
		t_error("ERROR: ATTEMPTING TO READ A CLOSED FILE");
	}
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
	
	