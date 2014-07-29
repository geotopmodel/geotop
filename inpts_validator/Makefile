#GEOtop inpts_validator
#Copyright (C) 2014 Exact Lab. srl

all: inpts_validator

definitions.h:
	cp definitions.h.in definitions.h

y.tab.c y.tab.h: geotop.y
	yacc -v -d geotop.y

geotop.lex.c: y.tab.h geotop.l definitions.h
	flex -o $@ geotop.l

inpts_validator: y.tab.c geotop.lex.c
	cc -o $@ y.tab.c geotop.lex.c

test: test_lex test_parser
	cat geotop.inpts | ./test_lex
	cat geotop.inpts | ./test_parser

test_parser: y.tab.c geotop.lex.c
	cc -DVERBOSE -o $@ y.tab.c geotop.lex.c

test_lex: geotop.lex.c
	cc -DDEBUG -lfl -o $@ geotop.lex.c

clean:
	-rm inpts_validator
	-rm y.output
	-rm test*
	-rm *c
	-rm *h

