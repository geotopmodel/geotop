#include "turtle.h"

STRINGBIN *read_plane_strings_restricted(FILE *inputfile,long length,long maxbuffersize);

STRINGBIN *read_plane_strings_restricted(FILE *inputfile,long length,long maxbuffersize)

{

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
			buffer=realloc(buffer,(buffer_size+1)*sizeof(char));
			if(buffer==NULL){
				t_error("Cannot expand the buffer");
			}
			buffer[buffer_size]='\0';
		} else {
			buffer_size+=BUFFERINCREMENT;
			buffer=realloc(buffer,(buffer_size+1)*sizeof(char));
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
