#include "turtle.h"
#include "t_datamanipulation.h"
#include "t_random.h"

#define MAX_NAME 256

/*--------------------------------------------------------------------------*/ 

/* Translating a buffer of longs into a list */

LONGPAIR* buffer2list(long  *s,LONGPAIR* nxt)


{
LONGPAIR* head;
if(s[0]==0 || s[1]==0)  
	return nxt; 
else {
	head=(LONGPAIR* )malloc(sizeof(LONGPAIR));
	if(!head) t_error("I cannot allocate further memory");
	head->i=s[0];
	head->j=s[1];
	head->next = buffer2list(s+2,nxt);
	return head;
	}
}


/*--------------------------------------------------------------------------*/ 

REALPAIR* realbuffer2list(double  *s,REALPAIR* nxt)


{

REALPAIR* head;
if(s[0]==0 || s[1]==0)  
	return nxt; 
else {
	head=(REALPAIR* )malloc(sizeof(REALPAIR));
	if(!head) t_error("I cannot allocate further memory");
	head->x=s[0];
	head->y=s[1];
	head->next = realbuffer2list(s+2,nxt);
	return head;
	}
}

/*--------------------------------------------------------------------------*/ 

LONGPAIR *new_longpair()
{

LONGPAIR *head;

	head=(LONGPAIR* )malloc(sizeof(LONGPAIR));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->i=0;
	head->j=0;

	return head;


}


/*--------------------------------------------------------------------------*/ 

LONGPOKER *new_longpoker()
{

LONGPOKER *head;

	head=(LONGPOKER* )malloc(sizeof(LONGPAIR));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->i=0;
	head->j=0;
	head->k=0;
	head->l=0;

	return head;


}

/*--------------------------------------------------------------------------*/ 

REALPAIR *new_realpair()
{

REALPAIR *head;

	head=(REALPAIR* )malloc(sizeof(REALPAIR));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->x=0;
	head->y=0;


	return head;


}
/*--------------------------------------------------------------------------*/ 

IX *new_ix()
{

IX *head;

	head=(IX* )malloc(sizeof(IX));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->i=0;
	head->x=0;


	return head;


}

/*--------------------------------------------------------------------------*/ 

IJX *new_ijx()
{

IJX *head;

	head=(IJX* )malloc(sizeof(IJX));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->i=0;
	head->j=0;
	head->x=0;


	return head;


}
/*--------------------------------------------------------------------------*/ 

XYZ *new_xyz()
{

XYZ *head;

	head=(XYZ* )malloc(sizeof(XYZ));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->x=0;
	head->y=0;
	head->z=0;


	return head;


}



PHRASE *new_word(long n,char *text)

{

long m;
PHRASE *head;
        if(n >  MAX_NAME) {
        	m= MAX_NAME;
        }else if(n <= 0){
            printf("\nWarning::returning a NULL PHRASE pointer\n");
        	return NULL;
        } else {
        	m=n;
        }
	head=(PHRASE*)malloc(sizeof(PHRASE));
	if(!head) t_error("I cannot allocate further memory");
	head->next = NULL;
	head->nh=m;
    head->word=(char *)malloc((m+1)*sizeof(char));
    strcpy(head->word, text);
        
	return head;


}

/*--------------------------------------------------------------------------*/ 


void print_longlist_elements(LONGPAIR * list,short columns)
{

long count=1;
LONGPAIR *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	putchar('\n');
	printf("%ld %ld\t",tmp->i,tmp->j);
	while(tmp->next!=NULL){
    	tmp=tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%ld %ld\t",tmp->i,tmp->j);
	}
	putchar('\n');
}
}


/*--------------------------------------------------------------------------*/ 


void print_pokerlist_elements(LONGPOKER * list,short columns)
{

long count=1;
LONGPOKER *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	putchar('\n');
	printf("%ld %ld %ld %ld\t",tmp->i,tmp->j,tmp->k,tmp->l);
	while(tmp->next!=NULL){
    	tmp=tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%ld %ld %ld %ld\t",tmp->i,tmp->j,tmp->k,tmp->l);
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 


void write_indexed_pokerlist_elements(FILE *ostream,LONGPOKER * list,short columns)
{

long count=1;
LONGPOKER *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	putchar('\n');
	fprintf(ostream,"%ld %ld %ld %ld\t",tmp->i,tmp->j,tmp->k,tmp->l);
	while(tmp->next!=NULL){
    	tmp=tmp->next;
    	count++;
    	if(count%columns==0) {putchar('\n');}
		fprintf(ostream,"%ld %ld %ld %ld\t",tmp->i,tmp->j,tmp->k,tmp->l);
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 

void print_reallist_elements(REALPAIR * list,short columns)
{

long count=1;

REALPAIR *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	putchar('\n');
	printf("%f %f\t",tmp->x,tmp->y);
	while(tmp->next!=NULL){
    	tmp=tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%f %f\t",tmp->x,tmp->y);
	}
	putchar('\n');
}
}
/*--------------------------------------------------------------------------*/ 

void print_ix_elements(IX * list,short columns)
{

long count=1;
IX *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	putchar('\n');
	printf("%ld %f\t",tmp->i,tmp->x);
	while(tmp->next!=NULL){
    	tmp=tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%ld %f\t",tmp->i,tmp->x);
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 

void print_ijx_elements(IJX * list,short columns)
{

long count=1;
IJX *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	putchar('\n');
	tmp=list;
	printf("%ld %ld %f\n", tmp->i, tmp->j, tmp->x);
	while( tmp->next!=NULL){
    	 tmp= tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%ld %ld %f\n", tmp->i, tmp->j, tmp->x);
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 

void print_xyz_elements(XYZ * list,short columns)
{

long count=1;
XYZ *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	putchar('\n');
	tmp=list;
	printf("%f %f %f\n", tmp->x, tmp->y, tmp->z);
	while( tmp->next!=NULL){
    	tmp= tmp->next;
    	count++;
    	if(count%columns==0) putchar('\n');
		printf("%f %f %f\n", tmp->x, tmp->y,tmp->z);
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 

void print_phrase_elements(PHRASE * list,short columns)

{

long count=1;
PHRASE *tmp;

if(list==NULL){ 
	printf(" NULL LIST\n");
}
else{
	tmp=list;
	if(columns==1){
		printf("%s\n", tmp->word);
		while( tmp->next!=NULL){
		    	tmp= tmp->next;
		    	count++;
			    printf("%s\n", tmp->word);	
	}
	}else{
	 	printf("%s ", tmp->word);
		while( tmp->next!=NULL){
    			tmp= tmp->next;
    			count++;	
    			if(count%columns==0) putchar('\n');
		        printf("%s ", tmp->word);
			
		}	
	}
	putchar('\n');
}
}

/*--------------------------------------------------------------------------*/ 

long count_longpair_elements(LONGPAIR * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_longpair_elements(head->next));
}

/*--------------------------------------------------------------------------*/ 

long count_longpoker_elements(LONGPOKER * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_longpoker_elements(head->next));
}

/*--------------------------------------------------------------------------*/ 

long count_ix_elements(IX * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_ix_elements(head->next));
}

/*--------------------------------------------------------------------------*/ 

long count_realpair_elements(REALPAIR * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_realpair_elements(head->next));
}

/*--------------------------------------------------------------------------*/ 

long count_ijx_elements(IJX * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_ijx_elements(head->next));
}


/*--------------------------------------------------------------------------*/ 

long count_xyz_elements(XYZ * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_xyz_elements(head->next));
}

/*--------------------------------------------------------------------------*/ 

long count_phrase_elements(PHRASE * head)

{
	if(head==NULL)
		return 0;
	else
		return(1+count_phrase_elements(head->next));
}



/*--------------------------------------------------------------------------*/ 
/* Pointing to the antecedent to a numbered element of the list.
If NULL the first is being pointed. */
LONGPAIR * point2longpair(LONGPAIR * head,long point)

{
   /* 
    printf("#%d#\n",point);
    print_longlist_elements(head,10); 
  */
	if(head==NULL) {
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1) {
        return head;
	}else {
	    point2longpair(head->next,point-1);
    }
 return NULL; 

}


/*--------------------------------------------------------------------------*/ 
LONGPOKER * point2longpoker(LONGPOKER * head,long point)

{
	if(head==NULL) {
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1) {
        return head;
	}else {
	    point2longpoker(head->next,point-1);
    }
 return NULL; 
}

/*--------------------------------------------------------------------------*/ 
REALPAIR * point2realpair(REALPAIR * head,long point)

{
    
	if(head==NULL) {
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1) {
        return head;
	}else {
	    point2realpair(head->next,point-1);
    }
 return NULL; 

}

/*--------------------------------------------------------------------------*/ 
IX * point2ix(IX * head,long point)

{
    
	if(head==NULL) {
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1){ 
        return head;
	}else {
	    point2ix(head->next,point-1);
    }
 return NULL; 

}

/*--------------------------------------------------------------------------*/ 
IJX * point2ijx(IJX * head,long point)

{
    
	if(head==NULL){ 
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1) {
        return head;
	}else {
	    point2ijx(head->next,point-1);
	 }   
 return NULL; 

}
/*--------------------------------------------------------------------------*/ 
XYZ * point2measure(XYZ * head,long point)

{
    
	if(head==NULL){ 
		t_error("This element does not exist in the list");
	}else if (point<=0){
		return NULL; 
    }else if (point==1) {
        return head;
	}else {
	    point2measure(head->next,point-1);
	}
	
	return NULL;

}

/*--------------------------------------------------------------------------*/ 

/* Deleting one element in a list and returning the pointer
to the head */


LONGPAIR * delete_longpair(LONGPAIR * head,LONGPAIR * ante)
{
LONGPAIR * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}



LONGPOKER * delete_longpoker(LONGPOKER * head,LONGPOKER * ante)
{
LONGPOKER * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}

/* Deleting one element in a list appending the first part of the list
preceding the element to the rest  */


LONGPAIR * deleteshuffle_longpair(LONGPAIR * head,LONGPAIR * ante)
{
LONGPAIR * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=NULL;
 head=appendto(q->next,head);
 
 }
free(q);
return head;
}

/*--------------------------------------------------------------------------*/ 

REALPAIR * delete_realpair(REALPAIR * head,REALPAIR * ante)
{

REALPAIR * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}
/*--------------------------------------------------------------------------*/ 

IX* delete_ix(IX * head,IX * ante)
{

IX * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}

/*--------------------------------------------------------------------------*/ 

IJX * delete_ijx(IJX * head,IJX * ante)
{

IJX * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}

/*--------------------------------------------------------------------------*/ 

XYZ * delete_xyz(XYZ * head,XYZ * ante)
{

XYZ * q;
if(head==NULL) t_error("NULL cannot be deleted");
if(ante==NULL) {
 q=head;
 head=head->next; 
 }
 else{
 q=ante->next;
 ante->next=q->next;
 }
free(q);
return head;
}


/*--------------------------------------------------------------------------*/ 

/*Selecting one point randomly in a list. It returns: NULL if the
point chosen is the first of the list. The antecedent to the point in any different case. 
The reason is that you need to know the antecedent for deleting the point from the list.*/

/* commentata per  warning: integer overflow in expression

LONGPAIR * select_longpair_randomly(LONGPAIR * pointer)
 
{

long length,utmp;
LONGPAIR * linktmp;

if(pointer==NULL) t_error("NULL cannot be selected");
length=count_longpair_elements(pointer);
utmp=rrand(length);
if(utmp>=length) t_error("wrong random number");
linktmp=point2longpair(pointer,utmp);
return linktmp;
}
/*--------------------------------------------------------------------------*/ 
/* commentata per  warning: integer overflow in expression
REALPAIR * select_realpair_randomly(REALPAIR * pointer)
 
{

unsigned long length,utmp;
REALPAIR * linktmp;

if(pointer==NULL) t_error("NULL cannot be selected");
length=count_realpair_elements(pointer);
utmp=rrand(length);
if(utmp>=length) t_error("wrong random number");
linktmp=point2realpair(pointer,utmp); 
return linktmp;
}

/*--------------------------------------------------------------------------*/ 

/* commentata per  warning: integer overflow in expression
IX * select_ix_randomly(IX * pointer)
 
{

unsigned long length,utmp;
IX * linktmp;

if(pointer==NULL) t_error("NULL cannot be selected");
length=count_ix_elements(pointer);
utmp=rrand(length);
if(utmp>=length) t_error("wrong random number");
linktmp=point2ix(pointer,utmp); 
return linktmp;
}


/*--------------------------------------------------------------------------*/ 

/* commentata per  warning: integer overflow in expression
IJX * select_ijx_randomly(IJX * pointer)
 
{

unsigned long length,utmp;
IJX * linktmp;

if(pointer==NULL) t_error("NULL cannot be selected");
length=count_ijx_elements(pointer);
utmp=rrand(length);
if(utmp>=length) t_error("wrong random number");
linktmp=point2ijx(pointer,utmp); 
return linktmp;
}

/*--------------------------------------------------------------------------*/ 
/* commentata per  warning: integer overflow in expression
XYZ * select_xyz_randomly(XYZ * pointer)
 
{

unsigned long length,utmp;
XYZ * linktmp;

if(pointer==NULL) t_error("NULL cannot be selected");
length=count_xyz_elements(pointer);
utmp=rrand(length);
if(utmp>=length) t_error("wrong random number");
linktmp=point2measure(pointer,utmp); 
return linktmp;

}


/* Recursive deletion of a list */
/*--------------------------------------------------------------------------*/ 

void delete_longpair_list(LONGPAIR * head)
{
if(head!=NULL){
	delete_longpair_list(head->next);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

void delete_longpoker_list(LONGPOKER * head)
{
if(head!=NULL){
	delete_longpoker_list(head->next);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

void delete_realpair_list(REALPAIR * head)
{
if(head!=NULL){
	delete_realpair_list(head->next);
	free(head);
}
}
/*--------------------------------------------------------------------------*/ 

void delete_ix_list(IX * head)
{
if(head!=NULL){
	delete_ix_list(head->next);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

void delete_ijx_list(IJX * head)

{

if(head!=NULL){
	delete_ijx_list(head->next);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

void delete_xyz_list(XYZ * head)
{
if(head!=NULL){
	delete_xyz_list(head->next);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

void delete_phrase(PHRASE *head)
{
/* char ch; */
if(head!=NULL){
	delete_phrase(head->next);
	free(head->word);
	free(head);
}
}

/*--------------------------------------------------------------------------*/ 

LONGPAIR * appendto(LONGPAIR * head, LONGPAIR * element)
{

LONGPAIR * tmp=head;

if(head==NULL) return element;
else{
	while(tmp->next!=NULL) tmp=tmp->next;
	tmp->next=element;
	return head;
}
}


/*--------------------------------------------------------------------------*/ 

PHRASE *  join_words(PHRASE * head, PHRASE* element)

{

PHRASE * tmp=head;

if(head==NULL) {
	return element;
}
else{
	while(tmp->next!=NULL) tmp=tmp->next;
	tmp->next=element;
	return head;
}
}

/*--------------------------------------------------------------------------*/ 

LONGPAIR * prependto(LONGPAIR * head, LONGPAIR * element)
{

LONGPAIR* home;

home=element;
if(head==NULL){ 	   
	return element;
}else if(element==NULL){
	return head;	
}else {
while(home->next!=NULL){
	printf("k\n");
	home=home->next;
}
home->next=head;
return home;
}

}
/*--------------------------------------------------------------------------*/ 
LONGPAIR *  rotateleft(LONGPAIR * head)

{

LONGPAIR * tmp=head;

if(head==NULL) return NULL;
else{
while(tmp->next!=NULL) tmp=tmp->next; 

tmp->next=head;
tmp=head->next;
head->next=NULL;
return tmp;
}
}

/*--------------------------------------------------------------------------
This rotate the list n times
---------------------------------------------------------------------------*/ 
LONGPAIR *  rotate(LONGPAIR * head,int n)
{

LONGPAIR * tmp=head,*back=NULL;
int m=0;

if(n==0) {
	return tmp;
}else if(head==NULL){ 
	return NULL;
}else{
	while(tmp->next!=NULL && m<n) {back=tmp;tmp=tmp->next; m++;} 
	if(m==n) {
		back->next=NULL;
		back=tmp;
		while(back->next!=NULL) back=back->next;
		back->next=head;
		return tmp;
		}
	else if(tmp->next==NULL)
		return rotate(head,n-m);
}

return NULL;

}

