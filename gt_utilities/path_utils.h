
/* GEOTRIVIALUtilities is an under-construction set of  C functions which supports i/o interface
   and other utilities for models like GEOtop
GEOTRIVIALUtilities Version 1.0

file geo_trivial_symbols.h

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
 * @file path_utils.h
 * @Author Gianfranco Gallizia (skyglobe83@gmail.com) 
 * @brief Utility functions for path handling
 */

#ifndef PATH_UTILS_H
#define PATH_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#if __STDC_VERSION__ >= 199901L
    const size_t MAX_PATH_LENGTH = 254;
#else
#define MAX_PATH_LENGTH 254
#endif

/**
 * @brief Checks if file exists
 * @param[in] filename the path of the file
 * @return 0 if filename doesn't exists, 1 if filename exists, 2 if
 *         filename is a directory, 3 if filename is a special file
 *         and -1 on error.
 */
int gt_fileExists(const char* filename);

/**
 * @brief Attempts to create a directory
 * @param[in] path the path of the directory
 * @return 1 if successful 0 on error.
 */
int gt_makeDirectory(const char* path);

/**
 * @brief removes the last element of a path
 *
 * E.G.: "/home/user/file" becomes "/home/user"
 *
 * @param[in] path the path to pop
 * @return a newly created string that must be freed with free() by the
 *         caller or NULL in case of error.
 */
char* gt_popPath(const char* path);

#ifdef __cplusplus
}
#endif
#endif
