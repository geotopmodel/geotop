/*
 * GEOtop inpts_validator
 * Copyright (C) 2014 eXact lab srl
 */

%{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "definitions.h"
#include "y.tab.h"

void printSymbol(char *s) {
#ifdef DEBUG
    printf(" %s", s);
#endif
}

#ifdef DEBUG
YYSTYPE yylval;
#endif

/*strdup test*/
#ifndef _BSD_SOURCE
char* strdup(const char* s)
{
    size_t len = strlen(s) + 1;
    char* output = calloc(size_of(char), len);

    if (output)
    {
        strcpy(output, s);
    }

    return output;

}
#endif

%}

%option noyywrap nodefault yylineno
%%

\n { printSymbol("EOL\n"); return EOL; }
[#!].* { printSymbol("COMMENT"); return COMMENT; }
"INS" { printSymbol("INTEGRATION"); yylval.i = INS; return INTEGRATION; }
"CUM" { printSymbol("INTEGRATION"); yylval.i = CUM; return INTEGRATION; }
"AVG" { printSymbol("INTEGRATION"); yylval.i = AVG; return INTEGRATION; }
[Tt][Rr][Uu][Ee] { printSymbol("BOOL"); yylval.i = TRUE; return BOOL; }
[Ff][Aa][Ll][Ss][Ee] { printSymbol("BOOL"); yylval.i = FALSE; return BOOL; }
[a-zA-Z][a-zA-Z0-9_]* { printSymbol("IDENT"); yylval.s = strdup(yytext); return IDENT; }
"[" { printSymbol("BRACKET_OP"); return BRACKET_OP; }
"]" { printSymbol("BRACKET_CL"); return BRACKET_CL; }
(1Ds|1Dp|2D|3D) { printSymbol("DIMENSION"); return DIMENSION; }
"=" { printSymbol("EQUAL_SIGN"); return EQUAL_SIGN; }
"/" { printSymbol("SLASH_SEP"); return SLASH_SEP; }
":" { printSymbol("COLON_SEP"); return COLON_SEP; }
"::" { printSymbol("O_KEY_SEP"); return O_KEY_SEP; }
"," { printSymbol("ARRAY_SEP"); return ARRAY_SEP; }
[0-9][0-9]*\.?[0-9]*([eE][-+]?[0-9]+)? { printSymbol("NUM"); yylval.d = atof(yytext); return NUM; }
\"[^\"\n]*\" { printSymbol("STRING"); yylval.s = strdup(yytext); return STRING; }
. { /* Do Nothing */ }

%%

