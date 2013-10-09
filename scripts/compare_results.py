#!/usr/bin/env python
# -*- coding: utf-8 -*-#
# @(#)compare_results.py
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
"""
Script to compare results from two different output directories or files.
"""

__author__ = 'Antonio Messina <antonio.s.messina@gmail.com>'
__docformat__ = 'reStructuredText'

import argparse
import filecmp
import os
import pandas

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

quiet = False
verbose = False

class NotATabFile(Exception):
    pass

def p_verbose(text, newline=True):
    if verbose:
        print text,
    if newline: print

def p_quiet(text, newline=True):
    if not quiet:
        print text,
    if newline: print

def indent(text, tabs=6):
    lines = [" "*tabs + i for i in str(text).split('\n')]
    return str.join('\n', lines)

def compare_tab_files(f1, f2):
    """
    Compare the data in two files and check if they are *almost*
    equal, depending on a delta which depends on the field name.
    """
    try:
        data_ok = pandas.read_csv(f1)
        data_new = pandas.read_csv(f2)
    except:
        print "Error reading files."
        return False
    if data_ok.shape[1] == 1:
        raise NotATabFile

    # check fields one by one
    equals = True
    for field in data_ok:
        p_quiet("  field '%s': " % field, newline=False)
        if data_ok[field].dtype == pandas.np.dtype('O'):
            if (data_ok[field] == data_new[field]).all():
                p_quiet(" OK")
            else:
                p_quiet(" FAIL")
                equals = False
        else:
            diff = (data_ok[field] - data_new[field]).apply(abs)
            if (diff == 0).all().all():
                p_quiet(" OK")
            else:
                p_quiet(" FAIL")
                p_verbose(indent(diff.describe()))
                equals = False
    return equals

def compare_map_files(f1, f2):
    data_ok = pandas.read_csv(f1, sep=' ', skiprows=6, header=None)
    data_new = pandas.read_csv(f2, sep=' ', skiprows=6, header=None)
    if data_ok.shape != data_new.shape:
        p_quiet("FAIL")
        return False
    try:
        diff = (data_ok - data_new).apply(abs)
    except (ValueError, TypeError):
        if not filecmp.cmp(f1, f2):
            p_quiet("FAIL")
            return False
    if (diff == 0).all().all():
        p_quiet("OK")
        return True
    else:
        p_quiet("FAIL")
        p_verbose(indent(diff.describe(), tabs=4))
        return False

def compare_files(f1, f2):
    p_quiet("Comparing '%s' and '%s':" % (f1, f2), newline=False)
    try:
        return compare_tab_files(f1, f2)
    except NotATabFile:
        return compare_map_files(f1, f2)

def build_file_pairs(path1, path2):
    """
    @returns (pairs, missing_in_1, missing_in_2)
    """
    files1 = []
    files2 = []
    if os.path.isdir(path1):
        if not os.path.isdir(path2):
            raise RuntimeError("cannot mix files and directories")
        files1 = os.listdir(path1)
        files2 = os.listdir(path2)
        common = [i for i in files1 if i in files2]
        missing1 = [i for i in files2 if i not in files1]
        missing2 = [i for i in files1 if i not in files2]
        files1 = [os.path.join(path1, i) for i in common]
        files2 = [os.path.join(path2, i) for i in common]
        return (zip(files1, files2),  missing1, missing2)
    else:
        return (((path1, path2),), (), ())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('path1')
    parser.add_argument('path2')
    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-q', '--quiet', action="store_true")
    args = parser.parse_args()
    # global quiet
    quiet = args.quiet
    # global verbose
    verbose = args.verbose
    if not os.path.exists(args.path1):
        raise RuntimeError("Invalid path: %s", args.path1)
    if not os.path.exists(args.path2):
        raise RuntimeError("Invalid path: %s", args.path2)

    (pairs, missing1, missing2) = build_file_pairs(args.path1, args.path2)

    for path in missing1:
        p_quiet("MISSING file %s/%s" % (args.path1, path))
    for path in missing2:
        p_quiet("MISSING file %s/%s" % (args.path2, path))

    results = []
    nres = len(pairs)
    it = 0
    for (f1, f2) in pairs:
        results.append(compare_files(f1, f2))
    print "Different files: %d" % results.count(False)
    print "Correct files: %d" % results.count(True)
    print "Missing files in %s: %d" % (args.path1, len(missing1))
    print "Missing files in %s: %d" % (args.path2, len(missing2))
