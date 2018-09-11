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
    if((fp=fopen(name,"r"))!=nullptr ){
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
    // free(basedir);
    exit(1);
  }

/*  free(basedir); */


  fp = fopen(name, mode);
  if (fp == nullptr) {
    printf("%s\n", name);
    t_error(" The specified file could not be opened ");

    return nullptr;

  } else {

    return fp;

  }

}


/**-----------------------------------------------------------------------*/
FILE *t_fclose(FILE *stream)
/* Safe file close */
{


  if (stream == nullptr) {
    printf(" An attemp was made to close an already closed file ");
  } else {
    fclose(stream);
  }

  return nullptr;

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









