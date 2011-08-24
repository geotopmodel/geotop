#include "turtle.h"

/*-----------------------------------------------------------------------*/
void t_error(char *error_text)
/* Error handling */
{

	/* void exit(); */
	char c;
	fprintf(stderr,"\nError::Run time error\n");
	fprintf(stderr,"Error::%s\n",error_text);
	fprintf(stderr,"........exiting...\n");
	c = getchar();
	exit(1);
}

