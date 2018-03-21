
/* GEOTRIVIALUtilities is an under-construction set of  C functions which
supports i/o interface and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

Copyright (c), 2011 Emanuele Cordano

This file is part of GEOTRIVIALUtilities.
 GEOTRIVIALUtilities is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GEOTRIVIALUtilities is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License v. 3.0 for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file path_utils.c
 * @Author Gianfranco Gallizia (skyglobe83@gmail.com)
 * @brief Utility functions for path handling
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "path_utils.h"

int gt_fileExists(const char *filename)
{
  int ret = 0, saved_errno;
  struct stat file_stat;

  ret = stat(filename, &file_stat);

  if (ret == 0)
    {
      // file exists
      if (S_ISREG(file_stat.st_mode))
        ret = 1;
      else if (S_ISDIR(file_stat.st_mode))
        ret = 2;
      else
        ret = 3;
    }
  else
    {
      // file doesn't exists or stat returned an error
      saved_errno = errno;
      if (saved_errno == 2)  // No such file or directory
        ret = 0;
      else
        ret = -1;
#ifdef DEBUG
      fprintf(stderr, "Error code: %d\nMessage: %s\n", saved_errno,
              strerror(saved_errno));
#endif
    }

  return ret;
}

char *gt_popPath(const char *path)
{
  char *ret = NULL;
  size_t len;
  int i, found = 0;

  len = strlen(path);

  if (len > MAX_PATH_LENGTH || len == 0) return NULL;

  for (i = len; i >= 0; --i)
    {
      if (path[i] == '/')
        {
          found = 1;
          break;
        }
    }

  if (!found) return NULL;

  if (i == 0) i++;  // Handles "/something" case

  ret = calloc(i + 1, sizeof(char));

  if (ret != NULL) strncpy(ret, path, i);

  return ret;
}

static int gt_recursiveMakeDirectory(char *path)
{
  int success = 0, ret = 1, saved_errno;
  size_t len;

  len = strlen((const char *)path);
  if (len > MAX_PATH_LENGTH || len == 0) return 0;

  ret = mkdir((const char *)path, 0755);
  saved_errno = errno;

  if (ret == 0) return 1;

  if (ret == -1 && saved_errno == ENOENT)
    {
      char *newpath = gt_popPath(path);
      if (newpath != NULL)
        {
          success = gt_recursiveMakeDirectory(newpath);
          free(newpath);
        }

      if (success)
        {
          ret = mkdir((const char *)path, 0755);

          if (ret == 0) success = 1;
        }
    }

  return success;
}

int gt_makeDirectory(const char *path)
{
  int ret = 0;
  char *newpath;

  switch (gt_fileExists(path))
    {
    case 0:
      // the directory doesn't exist
      // attempt to create it
      newpath = calloc(MAX_PATH_LENGTH + 1, sizeof(char));
      if (newpath != NULL)
        {
          strncpy(newpath, path, MAX_PATH_LENGTH);
          ret = gt_recursiveMakeDirectory(newpath);
          free(newpath);
        }
      else
        ret = 0;
      break;
    case 2:
      // file exists and is a directory
      // do nothing
      ret = 1;
      break;
    default:
      ret = 0;
      break;
    }

  return ret;
}

/*-----------------------------------------------------------------------*/
// Do we really need this function ?? TODO: remove it 15.12.2016

int mkdirp(const char *pathname, mode_t mode)
{
  /* From
   * http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/ */
  char parent[256], *p;

  /* make a parent directory path */
  strncpy(parent, pathname, sizeof(parent));
  parent[sizeof(parent) - 1] = '\0';
  for (p = parent + strlen(parent); *p != '/' && p != parent; p--)
    ;
  *p = '\0';

  /* try make parent directory */
  if (p != parent && mkdirp(parent, mode) != 0)
    {
#if defined(__CYGWIN__)
      return 0;
#else
      return -1;
#endif
    }

  /* make this one if parent has been made */
  if (mkdir(pathname, mode) == 0) return 0;

  /* if it already exists that is fine */
  if (errno == EEXIST)
    {
      struct stat sb;
      stat(pathname, &sb);
      if (!S_ISDIR(sb.st_mode))
        {
          /* pathname is NOT a directory! */
          return -1;
        }
      return 0;
    }
  return -1;
}
//----------------------------------

#ifdef __cplusplus
}
#endif
