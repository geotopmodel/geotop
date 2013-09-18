#!/usr/bin/env python
# -*- coding: utf-8 -*-#
# @(#)test_runner.py
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

__author__ = 'Antonio Messina <antonio.s.messina@gmail.com>'
__docformat__ = 'reStructuredText'

import filecmp
import os
import subprocess
import shutil
import sys
import tempfile

import unittest

TESTDIR = os.path.abspath(os.path.dirname(__file__))
GEOTOP = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../GEOtop_4'))

TESTS = ['example', 'small_example']

class TestValidRun(object):

    def setUp(self):
        os.chdir(TESTDIR)

    def compare_files(self, fpath_ok, fpath_new):
        """
        Compare the *content* of two files
        """
        assert filecmp.cmp(fpath_ok, fpath_new)

    def _test_template(self, directory):
        # Run geotop
        # import nose.tools; nose.tools.set_trace()
        os.chdir(directory)
        output = subprocess.check_output([GEOTOP, '.'], stderr=subprocess.STDOUT)

        assert os.path.isfile(os.path.join(directory, '_SUCCESSFUL_RUN'))
        
        resultdir = os.path.join(directory, 'results')

        for outputdir in os.listdir(resultdir):
            for fname in os.listdir(os.path.join(resultdir, outputdir)):
                fpath_ok = os.path.join(resultdir, outputdir, fname)
                fpath_new = os.path.join(directory, outputdir, fname)
                assert os.path.isfile(fpath_new)
                self.compare_files(fpath_ok, fpath_new)

    def test_generator(self):
        for d in TESTS:
            path = os.path.join(TESTDIR, d)
            if os.path.isdir(path):
                yield self._test_template, path

if __name__ == "__main__":
    import nose
    nose.run()
