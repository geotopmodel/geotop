/*
 * GEOtop inpts_validator
 * Copyright (C) 2014 Exact Lab. srl
 */

%{
#include <stdio.h>

extern int yylineno;

#include "definitions.h"

%}

%union {
    char* s;
    int i;
    double d;
}

/* Tokens */

%token EOL
%token COMMENT
%token IDENT
%token BRACKET_OP
%token BRACKET_CL
%token DIMENSION
%token INTEGRATION
%token EQUAL_SIGN
%token O_KEY_SEP
%token NUM
%token STRING
%token BOOL
%token ARRAY_SEP
%token SLASH_SEP
%token COLON_SEP

%%

configfile : /* nothing */
           | row EOL configfile
           | EOL configfile
           ;

row : section
    | comment
    | parameter
    ;

section : BRACKET_OP IDENT BRACKET_CL { VERBOSE_PRINT(" Section(%s)\n", yylval.s); free(yylval.s); yylval.s = NULL; }
        ;

comment : COMMENT
        ;

parameter : key EQUAL_SIGN value { VERBOSE_PRINT("\n"); }
          | key EQUAL_SIGN array { VERBOSE_PRINT("\n"); }
          ;

key : standard_key { VERBOSE_PRINT(" %s ", yylval.s); free(yylval.s); yylval.s = NULL; }
    | output_key { VERBOSE_PRINT(" Extended Key(%s) ", yyval.s); free(yylval.s); yylval.s = NULL; }
    ;

standard_key : IDENT
             ;

output_key : IDENT O_KEY_SEP DIMENSION O_KEY_SEP INTEGRATION
           ;

value : num_value
      | STRING { VERBOSE_PRINT(" String(%s) ", yylval.s); free(yylval.s); yylval.s = NULL; }
      | date { VERBOSE_PRINT(" Date "); }
      ;

num_value : BOOL { VERBOSE_PRINT(" Bool(%d) ", yylval.i); }
          | NUM { VERBOSE_PRINT(" Num(%f) ", yylval.d); }
          ;

array_value : num_value ARRAY_SEP
            ;

array : array_value
      | array array_value
      | array num_value
      ;

date : NUM SLASH_SEP NUM SLASH_SEP NUM NUM COLON_SEP NUM
     ;

%%

int main(int argc, char *argv[])
{
    yyparse();
    return 0;
}

yyerror(char *s)
{
    fprintf(stderr, "Error on line %d: %s\n", yylineno, s);
    exit(1);
}
