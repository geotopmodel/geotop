#include "turtle.h"

/*-----------------------------------------------------------------------*/
void t_error(const char *error_text)
/* Error handling */
{

  /* void exit(); */
  fprintf(stderr,"\nError::Run time error\n");
  fprintf(stderr,"Error::%s\n",error_text);
  fprintf(stderr,"........exiting...\n");
  exit(1);
}

