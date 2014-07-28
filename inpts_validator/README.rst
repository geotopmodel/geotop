GEOtop input file validator
===========================

This tool reads from standard input the contents of a geotop.inpts file
and exits with 0 if the file is valid or 1 if it's not and prints the
line number where it founded an error.

Usage
=====

To parse a ``geotop.inpts`` file::

    ./inpts_validator < geotop.inpts

Dependecies
===========

inpts_validator has been built using the following softwares:

- Flex 2.5.35
- Bison 2.4.1

You will also need a C compiler (tested with gcc 4.4.7 an glibc 2.12)
and make (tested with Gnu Make 3.81) to process the Makefile.

It should build with yacc/Berkeley Yacc and a POSIX-compliant make and
other C compilers but it hasn't been tested (caveat emptor).

Build
=====

Issuing a simple make without parameters should do the trick. If you
want you can run the tests with::

    make test

At the moment there is no install target in the Makefile.

