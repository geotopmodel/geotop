#ifndef T_LIST_H
#define T_LIST_H
#include "turtle.h"
#include "t_random.h"
#define MAX_NAME 256
/* LONGPAIR
Description:  It is a structure containing two long and a pointer to
a structure of the same type. */
struct  longpair{
	long i,j; 
	struct longpair *next;
	};
	
typedef struct longpair LONGPAIR;
/* LONGPOKER
Description:  It is a structure containing four long and a pointer to
a structure of the same type. 
*/
struct  longpoker{
	long i,j,k,l; 
	struct longpoker*next;
	};
	
typedef struct longpoker LONGPOKER;
/*  REAL PAIR
Description:  It is a structure containing two double and a pointer to
a structure of the same type. 
*/struct  realpair{
	double x,y; 
	struct realpair *next;
	};
	
typedef struct realpair REALPAIR;
/*  IX
Description:  It is a structure containing a long (i) and a double (x) and a pointer to
a structure of the same type.  
*/struct  ix{
	long i;
	double x; 
	struct ix *next;
	};
	
typedef struct ix IX;
/* IJX
Description:  It is a structure containing two long, double and a pointer to
a structure of the same type. 
 
*/struct  ijx{
	long i,j;
	double x; 
	struct ijx *next;
	};
	
typedef struct ijx IJX;
/* XYZ
Description:  It is a structure containing three double and a pointer to
a structure of the same type. 
 
*/struct  xyz{
	double x,y,z;
	struct xyz *next;
	};
	
typedef struct xyz XYZ;
struct measure{
 
long i,j;                     /*Pixel coordinates                                     */
DOUBLEVECTOR *data;           /*Pointing to a vector containing a set of measures     */
DOUBLEVECTOR *precision;      /*Pointing to a vector containing the precision of each 
                              measure. Usually a  vector gloabally defined         */
struct measure *next;};typedef struct measure MEASURES;
struct phrase{char  *word;
long   nh;
struct phrase *next; };
typedef struct phrase PHRASE;
/*Name: buffer2listVersion: 1.0Synopsis:LONGPAIR* buffer2list(long  *s,LONGPAIR* nxt);
REALPAIR* realbuffer2list(double  *s,REALPAIR* nxt);Authors & Date: Riccardo Rigon, 1998FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.cDescription:  They transform a vector of long or double (in the normal sense of
C programming) to a linear linked list of the correspondent type
 
Inputs a pointer to: 1) the vector to be transformed; 2) a LONG PAIR or REALPAIR
Return: the pointer to the head of the linear linked list created
See Also: Coupledfields_moments 
 
*/
LONGPAIR* buffer2list(long  *s,LONGPAIR* nxt);
REALPAIR* realbuffer2list(double  *s,REALPAIR* nxt);
/**Name: new_Version: 0.9Synopsis:LONGPAIR *new_longpair(void);
REALPAIR *new_realpair(void);
IX * new_ix(void);
IJX *new_ijx(void);
XYZ *new_xyz(void);
Description:  *_new_*creates aan element of a linear linked list of the specified types 
 Inputs: void
 
 Return: a pointer to: LONGPAIR,  REALPAIR,  IX and IJX respectively
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
LONGPAIR *new_longpair(void);
LONGPOKER *new_longpoker(void);
REALPAIR *new_realpair(void);
IX * new_ix(void);
IJX *new_ijx(void);
XYZ *new_xyz(void);
PHRASE *new_word(long n,char *text);
/**Name: print__elementsVersion: 1.0Synopsis:void print_longlist_elements(LONGPAIR * list,short columns);
void print_reallist_elements(REALPAIR * list,short columns);
void print_ix_elements(IX * list,short columns);
void print_ijx_elements(IJX * list,short columns);
void print_xyz_elements(XYZ * list,short columns);
Description:  They print the specified type of linear linked list
 
 Inputs: 1) the head of the linear linked list, the number of columns that
 will be printed on the screen
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
void print_pokerlist_elements(LONGPOKER * list,short columns);
void print_longlist_elements(LONGPAIR * list,short columns);
void print_reallist_elements(REALPAIR * list,short columns);
void print_ix_elements(IX * list,short columns);
void print_ijx_elements(IJX * list,short columns);
void print_xyz_elements(XYZ * list,short columns);
void print_phrase_elements(PHRASE * list,short columns);
/**Name: count__elementsVersion: 1.0Synopsis:long count_longpair_elements(LONGPAIR * head);
long count_realpair_elements(REALPAIR * head);
long count_ix_elements(IX * head);
long count_ijx_elements(IJX * head);
long count_xyz_elements(XYZ * head);
Description:  They counts the number of elements in the linear linked list of
the specified type
 
 Inputs: 1) The pointer to the head of the linear linked list
 Return: the number of elements in the list
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
long count_longpoker_elements(LONGPOKER * head);
long count_longpair_elements(LONGPAIR * head);
long count_realpair_elements(REALPAIR * head);
long count_ix_elements(IX * head);
long count_ijx_elements(IJX * head);
long count_xyz_elements(XYZ * head);
long count_phrase_elements(PHRASE * head);
void delete_phrase(PHRASE *head);
/**Name: delete_@Version: 1.0Synopsis:LONGPAIR * delete_longpair(LONGPAIR * head,LONGPAIR * ante);
REALPAIR * delete_realpair(REALPAIR * head,REALPAIR * ante);
IX * delete_ix(IX * head,IX * ante);
IJX * delete_ijx(IJX * head,IJX * ante);
XYZ * delete_xyz(XYZ * head,XYZ * ante);
Description:  They delete and deallocate an element of the  specified type from 
a linear linked list
 
 Inputs: a pointer to 1) The head of the linear linked list; 2) a pointer to
 the element preceding the element to be deleted (if this is Null the head is deleted)
 Return: the pointer to the head of the list
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
LONGPAIR * point2longpair(LONGPAIR * head,long point);
REALPAIR * point2realpair(REALPAIR * head,long point);
IX * point2ix(IX * head,long point);
IJX * point2ijx(IJX * head,long point);
XYZ * point2measure(XYZ * head,long point);
LONGPOKER * delete_longpoker(LONGPOKER * head,LONGPOKER * ante);
LONGPAIR * delete_longpair(LONGPAIR * head,LONGPAIR * ante);
REALPAIR * delete_realpair(REALPAIR * head,REALPAIR * ante);
IX * delete_ix(IX * head,IX * ante);
IJX * delete_ijx(IJX * head,IJX * ante);
XYZ * delete_xyz(XYZ * head,XYZ * ante);
LONGPOKER * point2longpoker(LONGPOKER * head,long point);LONGPAIR * deleteshuffle_longpair(LONGPAIR * head,LONGPAIR * ante);
/**Name: select_@_randomlyVersion: 1.0Synopsis:LONGPAIR * select_longpair_randomly(LONGPAIR * pointer);
REALPAIR * select_realpair_randomly(REALPAIR * pointer);
IX * select_ix_randomly(IX * pointer);
IJX * select_ijx_randomly(IJX * pointer);
XYZ * select_xyz_randomly(XYZ * pointer);
Description:  Select an element of a linear linked list choset at random
 
 Inputs: a pointer to the head of the list
 Return: the pointer to the selected element
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c  
*/
/* commentata per  warning: integer overflow in expression
LONGPAIR * select_longpair_randomly(LONGPAIR * pointer);
REALPAIR * select_realpair_randomly(REALPAIR * pointer);
IX * select_ix_randomly(IX * pointer);
IJX * select_ijx_randomly(IJX * pointer);
XYZ * select_xyz_randomly(XYZ * pointer);
*//**Name: delete__listVersion: 1.0Synopsis:void delete_longpair_list(LONGPAIR * head);
void delete_realpair_list(REALPAIR * head);
void delete_ix_list(IX * head);
void delete_ijx_list(IJX * head);
void delete_xyz_list(XYZ * head);
Description:  Delete a list of the specified type
 
 Inputs: A pointer to the head of the list
 Return: void
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
void delete_longpoker_list(LONGPOKER * head);
void delete_longpair_list(LONGPAIR * head);
void delete_realpair_list(REALPAIR * head);
void delete_ix_list(IX * head);
void delete_ijx_list(IJX * head);
void delete_xyz_list(XYZ * head);/**Name: appendtoVersion: 1.0Synopsis:LONGPAIR * appendto(LONGPAIR * head, LONGPAIR * element);Description:  appends element to a given linear linked list of the specified type
 
 Inputs: 1) a pointer to the head of the linear linked list; 2) a pointer to the element to be
 appended. 
 Return: a pointer to the head of the linear linked list
 
 Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
LONGPAIR * appendto(LONGPAIR * head, LONGPAIR * element);
PHRASE *  join_words(PHRASE * head, PHRASE* element);
/**Name: prependtoVersion: 1.0Synopsis:LONGPAIR * prependto(LONGPAIR * head, LONGPAIR * element);
Description:  append to the given element a the given linear linked list 
  Inputs: 1) a pointer to the head of the linear linked list; 2) a pointer to the element to be
 prepend 
 
 Return: a pointer to the head of the list
 Authors & Date: Riccardo Rigon, 1997
 
Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c  
*/
LONGPAIR * prependto(LONGPAIR * head, LONGPAIR * element);
/**Name: rotateleftVersion: 1.0Synopsis:LONGPAIR *  rotateleft(LONGPAIR * head);
Description:  Take a linear linked list and put the first element at the
end of the list
 Inputs: the pointer to the head of the list
 Return: a pointer to the new head of the list
 Authors & Date: Riccardo Rigon, 1997
 References: 
 Examples:  metropolis
 
Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
LONGPAIR *  rotateleft(LONGPAIR * head);
/**Name: rotateVersion: 1.0Synopsis:LONGPAIR *  rotate(LONGPAIR * head,int n);
Description:  appends the last element of a linear linked list to the head and
repeat the operation n times
 
 Inputs: the head of the linear linked list to be rotated and the number of rotations
 
 Return: the new head of the lsit
 
Authors & Date: Riccardo Rigon, 1997FILE: LIBRARIES/BASICS/t_list.h, LIBRARIES/BASICS/list.c
*/
LONGPAIR *  rotate(LONGPAIR * head,int n);
/** */
void write_indexed_pokerlist_elements(FILE *ostream,LONGPOKER * list,short columns);
#endif
