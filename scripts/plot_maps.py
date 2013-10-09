#!/usr/bin/env python
# -*- coding: utf-8 -*-#
# @(#)plot_maps.py
#
#
# Copyright (C) 2013, GC3, University of Zurich. All rights reserved.
#
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

__docformat__ = 'reStructuredText'

import os
import sys

import matplotlib.pylab as plt
import numpy as np

def read_file(fname):
    with open(fname) as fd:
        ncols = int(fd.readline().split()[-1])
        nrows = int(fd.readline().split()[-1])
        xllcorner = fd.readline().split()[-1]
        yllcorner = fd.readline().split()[-1]
        cellsize = fd.readline().split()[-1]
        NODATA_value = float(fd.readline().split()[-1])
        data = np.fromfile(fd, sep=' ')
    data.shape = (nrows, ncols)
    data = np.where(data == NODATA_value, np.nan, data)
    return data

def rename_file(fname, ext):
    """
    Function to rename of the ``fname`` file by appending an
    extension ``.N.old`` to the filename, and moving the ``ext``
    extension to the end of the path.

    For instance, create_backup("fname.png", ".png") should rename:
    "fname.png" => "fname.1.old.png"
    """
    if not fname.endswith(ext):
        print "WARNING: file %s does not end with %s" % (fname, ext)
        return None
    basename = fname[:-len(ext)]

    for i in range(1, 100):
        new_fname = "%s.%d%s" % (basename, i, ext)
        if os.path.exists(new_fname):
            continue
        os.rename(fname, new_fname)
        return new_fname
    return None

def plot_map_file(fname):
    data = read_file(fname)
    plt.clf()
    plt.imshow(data)
    plt.colorbar()
    pngfile = fname + '.png'
    if os.path.exists(pngfile):
        moved_to = rename_file(pngfile, '.png')
        if not moved_to:
            print "ERROR: file `%s` already exists and it was impossible to rename" % pngfile
            return False
    plt.savefig(pngfile)
    return True

if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv or len(sys.argv) == 1:
        prog = os.path.basename(sys.argv[0])
        print """
Usage: %s path [path2 [path3]]

This script accept as argument one or more path to a DEM file, and it
will produce a PNG picture for each file, named ``file.png``.

Example:

    %s dem.asc

will produce a picture ``dem.asc.png`` in the same directory as the
``dem.asc`` file.

If the destination file ``path.png`` already exists, it will be
renamed as ``path.N.png``, where N is a number between 1 and 99. If
more than 99 backup files exist, an error is printed.
""" % (prog, prog)
    for path in sys.argv[1:]:
        if os.path.isfile(path):
            print "Plotting file " + path + "...",
            if plot_map_file(path):
                print " done."
