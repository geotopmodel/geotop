
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
 * @file path_utils.c
 * @Author Gianfranco Gallizia (skyglobe83@gmail.com) 
 * @brief Utility functions for path handling
 */

#ifdef __cplusplus
extern "C" {
#endif

#include "path_utils.h"

int gt_fileExists(const char* filename)
{
    int ret = 0, saved_errno;
    struct stat file_stat;

    ret = stat(filename, &file_stat);

    if (ret == 0)
    {
        //file exists
        if (S_ISREG(file_stat.st_mode))
            ret = 1;
        else if (S_ISDIR(file_stat.st_mode))
            ret = 2;
        else 
            ret = 3;
    }
    else
    {
        //file doesn't exists or stat returned an error
        saved_errno = errno;
        if (saved_errno == 2)   //No such file or directory
            ret = 0;
        else
            ret = -1;
#ifdef DEBUG
        fprintf(stderr, "Error code: %d\nMessage: %s\n", saved_errno, strerror(saved_errno));
#endif
    }

    return ret;

}

#ifdef __cplusplus
}
#endif
