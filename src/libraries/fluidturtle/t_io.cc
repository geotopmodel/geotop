#include "turtle.h"
#include "t_utilities.h"

#include <sys/stat.h>
#include <libgen.h>
#include <limits.h>

extern const char *WORKING_DIRECTORY ;


int mkdirp(const char *pathname, mode_t mode) {
  /* From http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/ */
  char parent[PATH_MAX], *p;

  /* make a parent directory path */
  strncpy(parent, pathname, sizeof(parent));
  parent[sizeof(parent) - 1] = '\0';
  for (p = parent + strlen(parent); *p != '/' && p != parent; p--);
  *p = '\0';

  /* try make parent directory */
  if (p != parent && mkdirp(parent, mode) != 0)
    return -1;

  /* make this one if parent has been made */
  if (mkdir(pathname, mode) == 0)
    return 0;

  /* if it already exists that is fine */
  if (errno == EEXIST) {
    struct stat sb;
    stat(pathname, &sb);
    if (!S_ISDIR(sb.st_mode)) {
      /* pathname is NOT a directory! */
      return -1;
    }
    return 0;
  }
  return -1;
}

/**-----------------------------------------------------------------------*/
/* Safe file open */
FILE *t_fopen(const char *name, const char *mode) {

  //char newname[256];
  FILE *fp = nullptr;
  char *basedir = dirname(strdup(name));
  int ret = 0;

  /*if(strcmp(mode,"w")==0 || strcmp(mode,"wb")==0){
    if((fp=fopen(name,"r"))!=NULL ){
      // The file already exist
      //printf("\nWarning::Overwriting the file %s \n",name);
      strcpy(newname,name);
      strcat(newname,".old");
      t_fclose(fp);
      //rename(name,newname);
    }
  }*/

  ret = mkdirp(basedir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (-1 == ret) {
    fprintf(stderr, "ERROR: Unable to create parent directory `%s`. Exiting.\n",
            basedir);
    free(basedir);
    exit(1);
  }
/*  free(basedir); */

  fp = fopen(name, mode);
  if (fp == NULL) {
    printf("%s\n", name);
    t_error(" The specified file could not be opened ");

    return NULL;

  } else {

    return fp;

  }

}


/**-----------------------------------------------------------------------*/
FILE *t_fclose(FILE *stream)
/* Safe file close */
{


  if (stream == NULL) {
    printf(" An attemp was made to close an already closed file ");
  } else {
    fclose(stream);
  }

  return NULL;

}


/**-----------------------------------------------------------------------*/
long write_shortmatrix_elements(FILE *output, SHORTMATRIX *m, long maxcols)
/* Write a matrix of floating point numbers to a file */

{

  long tmp = 0, i, j;
  //long count=0;


  if (output == NULL) {
    t_error("The input file was not opened properly");
  } else if (m == NULL || m->co == NULL || m->isdynamic != 1) {
    t_error("The matrix was not allocated properly");
  } else if (m->nrl > m->nrh || m->ncl > m->nch) {
    t_error("The matrix has no proper dimensions");
  } else {

    for (i = m->nrl; i <= m->nrh; i++) {
      for (j = m->ncl; j <= m->nch; j++) {
        tmp = fprintf(output, "%hd\t ", m->co[i][j]);
        if (tmp == EOF) {
          printf("Error in storing data::Unespected End of file encountered\n");
          return EOF;
        }

        if (j % maxcols == 0 && j != m->nch) putc('\n', output);
      }
      putc('\n', output);
    }
    putchar('\n');
  }

  return OK;

}


/**-----------------------------------------------------------------------*/
long write_intmatrix_elements(FILE *output, INTMATRIX *m, long maxcols)
/* Write a matrix of int numbers to a file */

{

  long tmp = 0, i, j;
  //long count=0;

  if (output == NULL) {
    t_error("The input file was not opened properly");
  } else if (m == NULL || m->co == NULL || m->isdynamic != 1) {
    t_error("The matrix was not allocated properly");
  } else if (m->nrl > m->nrh || m->ncl > m->nch) {
    t_error("The matrix has no proper dimensions");
  } else {

    for (i = m->nrl; i <= m->nrh; i++) {
      for (j = m->ncl; j <= m->nch; j++) {
        tmp = fprintf(output, "%d ", m->co[i][j]);
        if (tmp == EOF) {
          printf("Error in storing data::Unespected End of file encountered\n");
          return EOF;
        }

        if (j % maxcols == 0 && j != m->nch) putc('\n', output);
      }
      putc('\n', output);
    }
    putchar('\n');
  }

  return OK;
}


/**-----------------------------------------------------------------------*/
long write_longmatrix_elements(FILE *output, LONGMATRIX *m, long maxcols)
/* Write a matrix of long point numbers to a file */

{

  long tmp = 0, i, j;
  //long count=0;


  if (output == NULL) {
    t_error("The input file was not opened properly");
  } else if (m == NULL || m->co == NULL || m->isdynamic != 1) {
    t_error("The matrix was not allocated properly");
  } else if (m->nrl > m->nrh || m->ncl > m->nch) {
    t_error("The matrix has no proper dimensions");
  } else {

    for (i = m->nrl; i <= m->nrh; i++) {
      for (j = m->ncl; j <= m->nch; j++) {
        /*printf("i(%d,%d)m(%d)\n",i,j,m->co[i][j]);*/
        tmp = fprintf(output, "%ld ", m->co[i][j]);
        if (tmp == EOF) {
          printf("Error in storing data::Unespected End of file encountered\n");
          return EOF;
        }

        if (j % maxcols == 0 && j != m->nch) putc('\n', output);
      }
      putc('\n', output);
    }
  }

  return OK;
}


/**-----------------------------------------------------------------------*/
long write_floatmatrix_elements(FILE *output, FLOATMATRIX *m, long maxcols)
/* Write a matrix of floating point numbers to a file */

{

  long tmp = 0, i, j;
  //long count=0;

  if (output == NULL) {
    t_error("The input file was not opened properly");
  } else if (m == NULL || m->co == NULL || m->isdynamic != 1) {
    t_error("The matrix was not allocated properly");
  } else if (m->nrl > m->nrh || m->ncl > m->nch) {
    t_error("The matrix has no proper dimensions");
  } else {

    for (i = m->nrl; i <= m->nrh; i++) {
      for (j = m->ncl; j <= m->nch; j++) {
        tmp = fprintf(output, "%f ", m->co[i][j]);
        if (tmp == EOF) {
          printf("Error in storing data::Unespected End of file encountered\n");
          return EOF;
        }

        if (j % maxcols == 0 && j != m->nch) putc('\n', output);
      }
      putc('\n', output);
    }
    putchar('\n');
  }

  return OK;

}


/**-----------------------------------------------------------------------*/
long write_doublematrix_elements(FILE *output, DOUBLEMATRIX *m, long maxcols)
/* Write a matrix of double numbers to a file */

{

  long tmp = 0, i, j;
  //long count=0;



  if (output == NULL) {
    t_error("The input file was not opened properly");
  } else if (m == NULL || m->co == NULL || m->isdynamic != 1) {
    t_error("The matrix was not allocated properly");
  } else if (m->nrl > m->nrh || m->ncl > m->nch) {
    t_error("The matrix has no proper dimensions");
  } else {

    for (i = m->nrl; i <= m->nrh; i++) {
      for (j = m->ncl; j <= m->nch; j++) {
        tmp = fprintf(output, "%g\t", m->co[i][j]);
        if (tmp == EOF) {
          printf("Error in storing data::Unespected End of file encountered\n");
          return EOF;
        }

        if (j % maxcols == 0 && j != m->ncl) putc('\n', output);
      }
      putc('\n', output);
    }
    /*putchar('\n');*/
  }

  return OK;

}


/**-------------------------------------------------------------------------*/

char *join_strings(const char *first, const char *second) {

  char *string;
  int len = strlen(first) + strlen(second) + 2;

  string = (char *) malloc(len * sizeof(char));
  string = strcpy(string, first);
  string = strcat(string, second);

  return string;
}



/**-----------------------------------------------------------------------*/









