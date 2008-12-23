
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
    
#include "turtle.h"
#include "import_ascii.h"
#include "extensions.h"

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

double *read_grassascii(double *header, double novalue, char *name, long max_figures)

{	
	FILE *f;
	long cont,i,j,ch[max_figures],r,c,nr,nc;
	short sgn, end=0;
	double *dtm;
	
	f=fopen(join_strings(name,ascii_grass),"r");
	if(f==NULL){
		printf("\nFile %s doesn't exist",join_strings(name,ascii_grass));
		t_error("Fatal error");
	}
	
	//read header
	for(i=0;i<=5;i++){	
		header[i]=0.0;
		cont=0;
		sgn=0;
		do{
			ch[0]=fgetc(f);
			error_message(1,ch[0],-1,10,0,name);
			cont+=1;
			if(cont==1){
				if( (i==0 && ch[0]!=110) || (i==1 && ch[0]!=115) || (i==2 && ch[0]!=101) || (i==3 && ch[0]!=119) || (i==4 && ch[0]!=114) || (i==5 && ch[0]!=99) ){
					printf("\n Warning: check if the file %s is in grass ascii format, the header is not ok \n",join_strings(name,ascii_grass));
					stop_execution();
				}
			}
		}while(ch[0]<=44 || ch[0]==47 || ch[0]>=58);
		cont=0;
		if(ch[0]==45){
			sgn=1;	//sign minus
			ch[0]=fgetc(f);
			if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
				printf("Error: - not followed by a number for r=%ld c=%ld, file %s",r,c,name,join_strings(name,ascii_grass));
				t_error("Fatal error");
			}
		}
		if(ch[0]>=48 && ch[0]<=57){		
			do{
				cont+=1;
				if(cont>=max_figures){
					printf("Increase max_figures_in_a_number in file %s",name,join_strings(name,ascii_grass));
					t_error("Fatal error");
				}
				ch[cont]=fgetc(f);
			}while(ch[cont]>=48 && ch[cont]<=57);
			for(j=0;j<=cont-1;j++){
				header[i]+=(ch[j]-48)*pow(10,cont-j-1);
			}
			ch[0]=ch[cont];
		}
		if(ch[0]==46){
			cont=0;
			do{
				ch[0]=fgetc(f);
				if(ch[0]>=48 && ch[0]<=57){
					cont+=1;
					header[i]+=(ch[0]-48)*pow(10,-cont);
				}
			}while(ch[0]>=48 && ch[0]<=57);
		}
		error_message(1,ch[0],-1,0,0,name);			
		if(ch[0]!=10){
			do{
				ch[0]=fgetc(f);
			}while(ch[0]!=10);
		}
		if(sgn==1) header[i]*=-1;
	}
	if(header[1]>=header[0]){
		printf("In file %s south larger than or equal to north",join_strings(name,ascii_grass));
		t_error("Fatal error");
	}
	if(header[3]>=header[2]){
		printf("In file %s west larger than or equal to east",join_strings(name,ascii_grass));
		t_error("Fatal error");
	}	
	if(header[4]<=0 || header[5]<=0){
		printf("In file %s nrows or ncols negative or null",join_strings(name,ascii_grass));
		t_error("Fatal error");
	}			
	
	//read matrix
	nr=(long)header[4];
	nc=(long)header[5];
	dtm = (double *) malloc(nr*nc*sizeof(double));
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			dtm[(r-1)*nc+c-1]=0.0;
		}
	}			
	for(r=1;r<=nr;r++){
		c=0;
		do{
			c+=1;
			sgn=0;
			do{
				ch[0]=fgetc(f);	
				if(ch[0]==10){	//end of line
					if(c==nc+1 || r>nr){
						r++;
						c=1;
					}else if(c!=1){
						printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_grass));
						t_error("Fatal error");
					}
				}
				if(ch[0]==-1){	//end of file
					if( ((r==nr)&&(c==nc+1)) || (r>nr) ){
						end=1;
						c=nc;
						r=nr;
					}else{
						printf("Number of rows less than declared in file %s",join_strings(name,ascii_grass));
						t_error("Fatal error");
					}
				}				
				if(ch[0]==58){
					printf("Header cannot consist of more than 6 lines, error in file %s",join_strings(name,ascii_grass));
					t_error("Fatal error");
				}
			}while(((ch[0]>=0)&&(ch[0]<=41))||(ch[0]==44)||(ch[0]==47)||(ch[0]>=58));
			//while ch is different from "*" "." "+" "-" or a number
			if(end==0){
				if(ch[0]==45) sgn=1;	//sign minus
				if(ch[0]==43 || ch[0]==45){
					ch[0]=fgetc(f);
					if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
						printf("Error: +- not followed by a number for r=%ld c=%ld, file %s",r,c,join_strings(name,ascii_grass));
						t_error("Fatal error");
					}
				}
				if(ch[0]>=48 && ch[0]<=57){
					cont=0;
					do{
						cont+=1;
						ch[cont]=fgetc(f);
					}while(ch[cont]>=48 && ch[cont]<=57);
					//while ch[cont] is a number
					for(j=0;j<=cont-1;j++){
						dtm[(r-1)*nc+c-1]+=(ch[j]-48)*pow(10,cont-j-1);
					}
					ch[0]=ch[cont];
				}
				if(ch[0]==46){
					cont=0;
					do{
						ch[0]=fgetc(f);
						if(ch[0]>=48 && ch[0]<=57){
							cont+=1;
							dtm[(r-1)*nc+c-1]+=(ch[0]-48)*pow(10,-cont);
						}
					}while(ch[0]>=48 && ch[0]<=57);
				}else if(ch[0]==42){
					dtm[(r-1)*nc+c-1]=novalue;
					ch[0]=fgetc(f);
				}else if(ch[0]!=32 && ch[0]!=10){
					ch[0]=fgetc(f);
				}
				if(sgn==1) dtm[(r-1)*nc+c-1]*=(-1);
			}
		}while(ch[0]!=10 && ch[0]!=-1);		
		
		if(c<nc && r<=nr){
			printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_grass));
			t_error("Fatal error");
		}
		if(c>nc && r<=nr){
			printf("Number of cols higher than declared in row %ld in file %s",r,join_strings(name,ascii_grass));
			t_error("Fatal error");
		}
		if(ch[0]==-1 && r<nr){
			printf("Number of rows less than declared in file %s",join_strings(name,ascii_grass));
			t_error("Fatal error");
		}
	}

	fclose(f);
	
	return(dtm);
		
}

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

double *read_esriascii(double *header, double novalue, char *name, long max_figures)

{	
	FILE *f;
	long cont,i,j,ch[max_figures],r,c,nr,nc;
	short sgn, end=0;
	double *dtm;
	
	f=fopen(join_strings(name,ascii_esri),"r");
	if(f==NULL){
		printf("\nFile %s doesn't exist",join_strings(name,ascii_esri));
		t_error("Fatal error");
	}
	
	//read header
	for(i=0;i<=4;i++){	
		header[i]=0.0;
		cont=0;
		sgn=0;
		do{
			ch[0]=fgetc(f);
			error_message(2,ch[0],-1,10,0,name);
			if(ch[0]==32 && cont==0){
			}else{
				cont+=1;
			}
			if(cont==1){
				if( (i==0 && ch[0]!=110) || (i==1 && ch[0]!=110) || (i==2 && ch[0]!=120) || (i==3 && ch[0]!=121) || (i==4 && ch[0]!=99) ){
					printf("\n Warning: check if the file %s is in esri ascii format, the header is not ok \n",join_strings(name,ascii_esri));
					stop_execution();
				}
			}			
			if( cont==5 && i==2 && ch[0]==101 ){
				printf("\n Warning: only xllcorner and yllcorner in the header of %s are allowed, if it is xllcenter and yllcenter the map cannot be correctly read \n",join_strings(name,ascii_esri));
				stop_execution();				
			}
		}while(ch[0]<=44 || ch[0]==47 || ch[0]>=58);
		cont=0;
		if(ch[0]==45){
			sgn=1;	//sign minus
			ch[0]=fgetc(f);
			if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
				printf("Error: +- not followed by a number for r=%ld c=%ld, file %s",r,c,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
		}
		if(ch[0]>=48 && ch[0]<=57){
			do{
				cont+=1;
				if(cont>=max_figures){
					printf("Increase max_figures_in_a_number in file %s",join_strings(name,ascii_esri));
					t_error("Fatal error");
				}
				ch[cont]=fgetc(f);			
			}while(ch[cont]>=48 && ch[cont]<=57);
			for(j=0;j<=cont-1;j++){
				header[i]+=(ch[j]-48)*pow(10,cont-j-1);
			}
			ch[0]=ch[cont];
		}
		if(ch[0]==46){
			cont=0;
			do{
				ch[0]=fgetc(f);
				if(ch[0]>=48 && ch[0]<=57){
					cont+=1;
					header[i]+=(ch[0]-48)*pow(10,-cont);
				}
			}while(ch[0]>=48 && ch[0]<=57);
		}
		error_message(2,ch[0],-1,0,0,name);			
		if(ch[0]!=10){
			do{
				ch[0]=fgetc(f);
			}while(ch[0]!=10);
		}
		if(sgn==1) header[i]*=-1;
	}
	if(header[0]<=0 || header[1]<=0){
		printf("In file %s nrows or ncols negative or null %ld %ld",join_strings(name,ascii_esri),header[0], header[1]);
		t_error("Fatal error");
	}		

	//read matrix
	nr=(long)header[1];
	nc=(long)header[0];
	dtm = (double *) malloc(nr*nc*sizeof(double));
	for(r=1;r<=nr;r++){
		for(c=1;c<=nc;c++){
			dtm[(r-1)*nc+c-1]=0.0;
		}
	}

	//check if a novalue line is present
	do{
		ch[0]=fgetc(f);	
		error_message(2,ch[0],-1,10,0,name);
	}while(ch[0]==32);
	if(ch[0]==78){	//character "N", NODATA_value
		header[5]=0.0;
		cont=0;
		sgn=0;
		do{
			ch[0]=fgetc(f);
			error_message(2,ch[0],-1,10,0,name);
		}while(ch[0]<=44 || ch[0]==47 || ch[0]>=58);
		if(ch[0]==45){
			sgn=1;	//sign minus
			ch[0]=fgetc(f);
			if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
				printf("Error: +- not followed by a number for r=%ld c=%ld, file %s",r,c,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
		}
		if(ch[0]>=48 && ch[0]<=57){		
			do{
				cont+=1;
				if(cont>=max_figures){
					printf("Increase max_figures_in_a_number in file %s",join_strings(name,ascii_esri));
					t_error("Fatal error");
				}
				ch[cont]=fgetc(f);
			}while(ch[cont]>=48 && ch[cont]<=57);
			for(j=0;j<=cont-1;j++){
				header[5]+=(ch[j]-48)*pow(10,cont-j-1);
			}
			ch[0]=ch[cont];
		}
		if(ch[0]==46){
			cont=0;
			do{
				ch[0]=fgetc(f);
				if(ch[0]>=48 && ch[0]<=57){
					cont+=1;
					header[5]+=(ch[0]-48)*pow(10,-cont);
				}
			}while(ch[0]>=48 && ch[0]<=57);
		}
		error_message(2,ch[0],-1,0,0,name);			
		if(ch[0]!=10){
			do{
				ch[0]=fgetc(f);
			}while(ch[0]!=10);
		}
		if(sgn==1) header[i]*=-1;
		
		for(r=1;r<=nr;r++){
			c=0;
			do{
				c+=1;
				sgn=0;
				do{
					ch[0]=fgetc(f);	
					if(ch[0]==10){	//end of line
						if(c==nc+1 || r>nr){
							r++;
							c=1;
						}else if(c!=1){
							printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}
					if(ch[0]==-1){	//end of file
						if( ((r==nr)&&(c==nc+1)) || (r>nr) ){
							end=1;
							c=nc;
							r=nr;
						}else{
							printf("Number of rows less than declared in file %s",join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}				
				}while(((ch[0]>=0)&&(ch[0]<=41))||(ch[0]==44)||(ch[0]==47)||(ch[0]>=58));
				//while ch is different from "*" "." "+" "-" or a number
				if(end==0){
					if(ch[0]==45) sgn=1;	//sign minus
					if(ch[0]==43 || ch[0]==45){
						ch[0]=fgetc(f);
						if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
							printf("Error: +- not followed by a number for r=%ld c=%ld, file %s",r,c,join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}
					if(ch[0]>=48 && ch[0]<=57){
						cont=0;
						do{
							cont+=1;
							ch[cont]=fgetc(f);
						}while(ch[cont]>=48 && ch[cont]<=57);
						//while ch[cont] is a number
						for(j=0;j<=cont-1;j++){
							dtm[(r-1)*nc+c-1]+=(ch[j]-48)*pow(10,cont-j-1);
						}
						ch[0]=ch[cont];
					}
					if(ch[0]==46){
						cont=0;
						do{
							ch[0]=fgetc(f);
							if(ch[0]>=48 && ch[0]<=57){
								cont+=1;
								dtm[(r-1)*nc+c-1]+=(ch[0]-48)*pow(10,-cont);
							}
						}while(ch[0]>=48 && ch[0]<=57);
					}else if(ch[0]==42){
						dtm[(r-1)*nc+c-1]=novalue;
						ch[0]=fgetc(f);
					}else if(ch[0]!=32 && ch[0]!=10){
						ch[0]=fgetc(f);
					}
					if(sgn==1) dtm[(r-1)*nc+c-1]*=(-1);
					if(dtm[(r-1)*nc+c-1]==header[5]) dtm[(r-1)*nc+c-1]=novalue;
				}
			}while(ch[0]!=10 && ch[0]!=-1);		
		
			if(c<nc && r<=nr){
				printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
			if(c>nc && r<=nr){
				printf("Number of cols higher than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
			if(ch[0]==-1 && r<nr){
				printf("Number of rows less than declared in file %s",join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
		}

	}else{
				
		for(r=1;r<=nr;r++){
			c=0;
			do{
				c+=1;
				do{
					if(c==1 && r==1){
						sgn=0;
					}else{
						ch[0]=fgetc(f);	
						sgn=0;
					}
					if(ch[0]==10){	//end of line
						if(c==nc+1 || r>nr){
							r++;
							c=1;
						}else if(c!=1){
							printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}
					if(ch[0]==-1){	//end of file
						if( ((r==nr)&&(c==nc+1)) || (r>nr) ){
							end=1;
							c=nc;
							r=nr;
						}else{
							printf("Number of rows less than declared in file %s",join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}				
				}while(((ch[0]>=0)&&(ch[0]<=41))||(ch[0]==44)||(ch[0]==47)||(ch[0]>=58));
				//while ch is different from "*" "." "+" "-" or a number
				if(end==0){
					if(ch[0]==45) sgn=1;	//sign minus
					if(ch[0]==43 || ch[0]==45){
						ch[0]=fgetc(f);
						if(ch[0]<=45 || ch[0]==47 || ch[0]>=58){
							printf("Error: +- not followed by a number for r=%ld c=%ld, file %s",r,c,join_strings(name,ascii_esri));
							t_error("Fatal error");
						}
					}
					if(ch[0]>=48 && ch[0]<=57){
						cont=0;
						do{
							cont+=1;
							ch[cont]=fgetc(f);
						}while(ch[cont]>=48 && ch[cont]<=57);
						//while ch[cont] is a number
						for(j=0;j<=cont-1;j++){
							dtm[(r-1)*nc+c-1]+=(ch[j]-48)*pow(10,cont-j-1);
						}
						ch[0]=ch[cont];
					}
					if(ch[0]==46){
						cont=0;
						do{
							ch[0]=fgetc(f);
							if(ch[0]>=48 && ch[0]<=57){
								cont+=1;
								dtm[(r-1)*nc+c-1]+=(ch[0]-48)*pow(10,-cont);
							}
						}while(ch[0]>=48 && ch[0]<=57);
					}else if(ch[0]==42){
						dtm[(r-1)*nc+c-1]=novalue;
						ch[0]=fgetc(f);
					}else if(ch[0]!=32 && ch[0]!=10){
						ch[0]=fgetc(f);
					}
					if(sgn==1) dtm[(r-1)*nc+c-1]*=(-1);
				}
			}while(ch[0]!=10 && ch[0]!=-1);		
		
			if(c<nc && r<=nr){
				printf("Number of cols less than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
			if(c>nc && r<=nr){
				printf("Number of cols higher than declared in row %ld in file %s",r,join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
			if(ch[0]==-1 && r<nr){
				printf("Number of rows less than declared in file %s",join_strings(name,ascii_esri));
				t_error("Fatal error");
			}
		}
	}

	fclose(f);
	
	return(dtm);
		
}

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
	
void error_message(short format, long n, long n1, long n2, long n3, char *name)
//format=1 grassascii
//format=2 esriascii

{
	if(n==n1 || n==n2 || n==n3){
		if(format==1) printf("File %s incompleted, end of file or end of line reached",join_strings(name,ascii_grass));
		if(format==2) printf("File %s incompleted, end of file or end of line reached",join_strings(name,ascii_esri));
		t_error("Fatal error");
	}
}

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
