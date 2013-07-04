#include "turtle.h"
#include "t_utilities.h"
#include "math.h"

 /*-----------------------------------------------------------------------*/



long search_named_array(FILE *inputfile,char* string)

{


char ch;
//const char ocurl='{', ccurl='}';
//short ind=1;
//long buffer_index=0,buffer_size=0,blocksnumber=0,
long position;
//extern t_keywords T_KEYWORDS;
//FLOATVECTOR *vec=NULL;
HEADER h;




if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
}


ch=getc(inputfile);
while(ch!=EOF){
	position=simplefind(inputfile,":");
	ch=getc(inputfile);
	while(!isspace(ch) && !iscntrl(ch)){
		position-=1;
		fseek(inputfile,position,SEEK_SET);
		ch=getc(inputfile);
	}
	header_scan(inputfile,&h);
	if(strcmp(h.name,string)==0){
		fseek(inputfile,position,SEEK_SET);
		return position;
	}
}
t_error("The searched array does not exists in file");
return -1;
}


 /*-----------------------------------------------------------------------*/



long search_array(FILE *inputfile,char* type,char* category)

{







char ch;

//const char ocurl='{', ccurl='}';

//short ind=1;

//long buffer_index=0,buffer_size=0,blocksnumber=0,
long position;
extern t_keywords T_KEYWORDS;


//FLOATVECTOR *vec=NULL;

HEADER h;




if(inputfile==NULL){
	t_error("You tried to read from a closed file ");
}

ch=getc(inputfile);

while(ch!=EOF){

	position=simplefind(inputfile,":");

    ch=getc(inputfile);
	while(!isspace(ch) && !iscntrl(ch)){
		position-=1;
		fseek(inputfile,position,SEEK_SET);
		ch=getc(inputfile);
	}
    
	header_scan(inputfile,&h);
	if(strcmp(T_KEYWORDS.type[h.type],type)==0 && strcmp(T_KEYWORDS.category[h.category],category)==0){
		fseek(inputfile,position,SEEK_SET);
		return position;
	}

}

t_error("The searched array does not exists in file");
return -1;


}


/*-------------------------------------------------------------------------*/
long simplefind(FILE *istream,const char *string)

{


long position=0,stringlength=0,counter=0;
char ch='a', *newstring;

stringlength=strlen(string);
newstring=(char *)malloc((stringlength+1)*sizeof(char));
if(newstring==NULL) t_error("Not enough memory to allocate this string");

newstring[stringlength]='\0';


while(ch!= EOF){

	ch=getc(istream);
	if(ch==EOF) return EOF;
	if(ch==string[0]) {
		position=ftell(istream);
	    newstring[0]=ch;
	    counter=1;
	    while(counter < stringlength){
			newstring[counter]=getc(istream);
			if(newstring[counter]==EOF) return EOF;
			counter++;
		}
   newstring[counter]='\0';
		if(strcmp(newstring,string)==0){
			break;
		}else {
		  position++;
		  fseek(istream,position,SEEK_SET);
		}
		
	}
}

if(ch==EOF){
	return EOF;
}else {

	fseek(istream,position-1,SEEK_SET);	
	return position--;	
}


}


/*--------------------------------------------------------------------------------------*/


void meter( long  index, long rows,  short frequence,const char* message, const char* separator)

{
short m;
if(frequence>rows){
	printf("%s %ld/%ld%s", message, index,rows,separator);
}else{
	m=ceil (rows/frequence);
	if( (index  % m )==0){
		printf("%s %ld/%ld%s", message, index,rows,separator);
	}
}
} 

/*-------------------------------------------------------------------------*/

long ssimplefindkeyword(char *buffer,long bufferlength, const char *string)

{


long position=0,stringlength=0,counter=0,pos=0,n;
//char ch='a';
char* newstring;

stringlength=strlen(string);
newstring=(char *)malloc((stringlength+1)*sizeof(char));
if(newstring==NULL) t_error("Not enough memory to allocate this string");
newstring[stringlength]='\0';
newstring[0]='a';
while( position <=bufferlength){
	counter=0;
	while(isspace(buffer[position] && position <= bufferlength)){
		position++;
	}
	newstring[counter]=toupper(buffer[position]);
	position++;
	if(newstring[counter]==string[0]) {
	    pos=position-1;
	    counter++;
	    while(counter < stringlength && position <= bufferlength){
			if(isspace(buffer[position] )){
				position++;
			}else{
				newstring[counter]=toupper(buffer[position]);
				position++;
				counter++;						
				if(newstring[counter-1]!=string[counter-1]) break;
			}
						
		}
   		newstring[counter]='\0';
		if(strcmp(newstring,string)==0){
			break;
		}else {
		  position++;
		}
		
	}
}

if( position > bufferlength){
	return -1;
}else {

    for(n=pos;n<position;n++){
        buffer[n]='\0';
    }
    n=pos-1;
    while(iscntrl(buffer[n]) || isspace(buffer[n])){
    	buffer[n]='\0';	
    	n--;
    }
    return position;	
}


}

/*-------------------------------------------------------------------------------*/
short file_copy(FILE *destination,FILE *origin)
{

char ch;
long pos=0;


if(destination==NULL || origin==NULL) t_error("This file was not opened yet");

		ch=getc(origin);
		while(ch!=EOF){
		       if(ch=='\n'){
		    		fprintf(destination,"\n");
		    		skip_whitespaces(origin); 
		    	}else{
		    	
		    		fprintf(destination,"%c",ch);
		    	}
			ch=getc(origin);
		}
pos=ftell(origin);
if(pos==0){
	printf("\nWarning::Empty file encountered");
	return 0;
} else {
	return 1;
}

	return 0;
}

/*-------------------------------------------------------------------------------*/
short join_strings_into(char *string, char *first,char *second)

{

	short len=0;
	
	
	len=strlen(first);
	len+=strlen(second)+2;
        if(len > FL)  t_error("Maximum directory length exceeded");
	string=strcpy(string,first);
	string=strcat(string,second);
	
	return 1;

}

/*-------------------------------------------------------------------------------*/
short join_3strings_into(char *string, char *first,char *second,char* third)

{

	short len=0;
	
	len=strlen(first);
	len+=strlen(second)+2;
	len+=strlen(third);
        if(len > FL)  t_error("Maximum directory length exceeded");
	string=strcpy(string,first);
	string=strcat(string,second);
	string=strcat(string,third);
	
	return 1;

}

/*--------------------------------------------------------------------------*/
void stop_execution(void)

{

char ch;

printf("\nPRESS RETURN TO CONTINUE\n");
scanf("%c",&ch);

}

/*--------------------------------------------------------------------------*/

void time2date(float time, long *giulian, long *year, long *month, long *day, long *hour, long *min, float *sec)

/* given a inputs 
	1:the time in second 
	2:date (giulian day, year, month, day, hour, min, sec)
   return as outputs the date updated for the time given
   	4:date (giulian day, year, month, day, hour, min, sec)
   bug: time have to be less than 1 year */

{
long dmon[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
long dyear,days,giu;
float t;

t=time;

if(fmod(*year,4)>0){
	dmon[2]=28;
	dyear=365;
}else{
	dmon[2]=29;
	dyear=366;
}

giu=*giulian+floor((time+*hour*3600+*min*60+*sec)/(24*3600));
/* giu=*giulian+floor((time)/(24*3600)); */
/* if changes the year the origin of time is shifted to
	the same hour of the January 1 of the next year*/
if(giu>dyear){ 
	t=time-(dyear-*giulian+1)*24*3600;
	*month=1;
	*day=1;
	*giulian=1;
	giu=giu-dyear;
	*year=*year+1;
	if(fmod(*year,4)>0){
		dmon[2]=28;
		dyear=365;
	}else{
		dmon[2]=29;
		dyear=366;
	}
}

days=0;
*month=0;
while(days<giu){
	*month=*month+1;	
	days+=dmon[*month];
}

days-=dmon[*month];
*day=giu-days;

/* shift the origin of time to the begin of the last day */
t=t-((((giu-*giulian-1)*24+(23-*hour))*60+(59-*min))*60+(60-*sec));
*hour=floor(t/3600);
*min=floor((t-*hour*3600)/60);
*sec=t-(*min+*hour*60)*60;
*giulian=giu;

}

/*--------------------------------------------------------------------------*/

void giulian2day(long giulian, long year, long *month, long *day)

/* given a inputs 
	giulian day, year
   return as outputs the date
   	month, day */

{
long dmon[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
long days;

if(fmod(year,4)>0){
	dmon[2]=28;
}else{
	dmon[2]=29;
}

days=0;
*month=0;
while(days<giulian){
	*month=*month+1;	
	days+=dmon[*month];
}

days-=dmon[*month];
*day=giulian-days;
}

/*--------------------------------------------------------------------------*/

void day2giulian(long year, long month, long day, long *giulian)

/* given a inputs 
	day, year, month
   return as outputs the julian day */

{
long dmon[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
long i;

if(fmod(year,4)>0){
	dmon[2]=28;
}else{
	dmon[2]=29;
}

*giulian=day;
for(i=1;i<month;i++){
	*giulian=*giulian+dmon[i];
}

}
	
/*--------------------------------------------------------------------------*/

double Fmin(double a, double b){
	
	double min=a;
	
	if(b<a) min=b;
	
	return(min);
	
}

/*--------------------------------------------------------------------------*/

double Fmax(double a, double b){

	double max=a;
	
	if(b>a) max=b;
	
	return(max);
	
}

/*--------------------------------------------------------------------------*/
