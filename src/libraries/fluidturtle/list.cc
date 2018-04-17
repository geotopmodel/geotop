#include "turtle.h"

#define MAX_NAME 256

/*--------------------------------------------------------------------------*/

/* Translating a buffer of longs into a list */

LONGPAIR *buffer2list(long  *s,LONGPAIR *nxt)


{
  LONGPAIR *head;
  if (s[0]==0 || s[1]==0)
    return nxt;
  else
    {
      head=(LONGPAIR * )malloc(sizeof(LONGPAIR));
      if (!head) t_error("I cannot allocate further memory");
      head->i=s[0];
      head->j=s[1];
      head->next = buffer2list(s+2,nxt);
      return head;
    }
}


/*--------------------------------------------------------------------------*/

REALPAIR *realbuffer2list(double  *s,REALPAIR *nxt)


{

  REALPAIR *head;
  if (s[0]==0 || s[1]==0)
    return nxt;
  else
    {
      head=(REALPAIR * )malloc(sizeof(REALPAIR));
      if (!head) t_error("I cannot allocate further memory");
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

  head=(LONGPAIR * )malloc(sizeof(LONGPAIR));
  if (!head) t_error("I cannot allocate further memory");
  head->next = NULL;
  head->i=0;
  head->j=0;

  return head;


}


/*--------------------------------------------------------------------------*/

long count_longpair_elements(LONGPAIR *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_longpair_elements(head->next));
}

/*--------------------------------------------------------------------------*/

long count_longpoker_elements(LONGPOKER *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_longpoker_elements(head->next));
}

/*--------------------------------------------------------------------------*/

long count_ix_elements(IX *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_ix_elements(head->next));
}

/*--------------------------------------------------------------------------*/

long count_realpair_elements(REALPAIR *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_realpair_elements(head->next));
}

/*--------------------------------------------------------------------------*/

long count_ijx_elements(IJX *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_ijx_elements(head->next));
}


/*--------------------------------------------------------------------------*/

long count_xyz_elements(XYZ *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_xyz_elements(head->next));
}

/*--------------------------------------------------------------------------*/

long count_phrase_elements(PHRASE *head)

{
  if (head==NULL)
    return 0;
  else
    return (1+count_phrase_elements(head->next));
}



/*--------------------------------------------------------------------------*/
/* Pointing to the antecedent to a numbered element of the list.
If NULL the first is being pointed. */
LONGPAIR *point2longpair(LONGPAIR *head,long point)

{
  /*
   printf("#%d#\n",point);
   print_longlist_elements(head,10);
  */
  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2longpair(head->next,point-1);
    }
  return NULL;

}


/*--------------------------------------------------------------------------*/
LONGPOKER *point2longpoker(LONGPOKER *head,long point)

{
  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2longpoker(head->next,point-1);
    }
  return NULL;
}

/*--------------------------------------------------------------------------*/
REALPAIR *point2realpair(REALPAIR *head,long point)

{

  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2realpair(head->next,point-1);
    }
  return NULL;

}

/*--------------------------------------------------------------------------*/
IX *point2ix(IX *head,long point)

{

  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2ix(head->next,point-1);
    }
  return NULL;

}

/*--------------------------------------------------------------------------*/
IJX *point2ijx(IJX *head,long point)

{

  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2ijx(head->next,point-1);
    }
  return NULL;

}
/*--------------------------------------------------------------------------*/
XYZ *point2measure(XYZ *head,long point)

{

  if (head==NULL)
    {
      t_error("This element does not exist in the list");
    }
  else if (point<=0)
    {
      return NULL;
    }
  else if (point==1)
    {
      return head;
    }
  else
    {
      point2measure(head->next,point-1);
    }

  return NULL;

}

/*--------------------------------------------------------------------------*/




/* Recursive deletion of a list */
/*--------------------------------------------------------------------------*/

void delete_longpair_list(LONGPAIR *head)
{
  if (head!=NULL)
    {
      delete_longpair_list(head->next);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

void delete_longpoker_list(LONGPOKER *head)
{
  if (head!=NULL)
    {
      delete_longpoker_list(head->next);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

void delete_realpair_list(REALPAIR *head)
{
  if (head!=NULL)
    {
      delete_realpair_list(head->next);
      free(head);
    }
}
/*--------------------------------------------------------------------------*/

void delete_ix_list(IX *head)
{
  if (head!=NULL)
    {
      delete_ix_list(head->next);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

void delete_ijx_list(IJX *head)

{

  if (head!=NULL)
    {
      delete_ijx_list(head->next);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

void delete_xyz_list(XYZ *head)
{
  if (head!=NULL)
    {
      delete_xyz_list(head->next);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

void delete_phrase(PHRASE *head)
{
  /* char ch; */
  if (head!=NULL)
    {
      delete_phrase(head->next);
      free(head->word);
      free(head);
    }
}

/*--------------------------------------------------------------------------*/

LONGPAIR *appendto(LONGPAIR *head, LONGPAIR *element)
{

  LONGPAIR *tmp=head;

  if (head==NULL) return element;
  else
    {
      while (tmp->next!=NULL) tmp=tmp->next;
      tmp->next=element;
      return head;
    }
}


/*--------------------------------------------------------------------------
This rotate the list n times
---------------------------------------------------------------------------*/
LONGPAIR   *rotate(LONGPAIR *head,int n)
{

  LONGPAIR *tmp=head,*back=NULL;
  int m=0;

  if (n==0)
    {
      return tmp;
    }
  else if (head==NULL)
    {
      return NULL;
    }
  else
    {
      while (tmp->next!=NULL && m<n) {back=tmp; tmp=tmp->next; m++;}
      if (m==n)
        {
          back->next=NULL;
          back=tmp;
          while (back->next!=NULL) back=back->next;
          back->next=head;
          return tmp;
        }
      else if (tmp->next==NULL)
        return rotate(head,n-m);
    }

  return NULL;

}

