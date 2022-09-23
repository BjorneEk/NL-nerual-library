/*==========================================================*
 *
 * @author Gustaf Franzén :: https://github.com/BjorneEk;
 *
 * assert and warning macros
 *
 *==========================================================*/

#ifndef _ASSERT_H_
#define _ASSERT_H_

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "types.h"
#include "macro.h"

#define log_err_none(err) {\
    fprintf(stderr, "[\033[31;1;4mError\033[0m]: %s | %s\n",(err), strerror( errno ));\
}
#define log_err_string(err, str) {\
    fprintf(stderr, "[\033[31;1;4mError\033[0m]: %s: '%s' | %s\n",(err), (str), strerror( errno ));\
}
#define log_err_int(err, i) {\
    fprintf(stderr, "[\033[31;1;4mError\033[0m]: %s: %i | %s\n",(err), (i), strerror( errno ));\
}
#define log_err_float(err, f) {\
    fprintf(stderr, "[\033[31;1;4mError\033[0m]: %s: %f | %s\n",(err), (f), strerror( errno ));\
}
#define log_err_char(err, c) {\
    fprintf(stderr, "[\033[31;1;4mError\033[0m]: %s: '%c' | %s\n",(err), (c), strerror( errno ));\
}
#define log_warn_none(warn) {\
  printf("[\033[33;1;4mWarning\033[0m]: %s\n", (warn));\
}
#define log_warn_string(warn, str) {\
  printf("[\033[33;1;4mWarning\033[0m]: %s: '%s'\n", (warn), (str));\
}
#define log_warn_int(warn, i) {\
  printf("[\033[33;1;4mWarning\033[0m]: %s: %i\n", (warn), (i));\
}
#define log_warn_float(warn, f) {\
  printf("[\033[33;1;4mWarning\033[0m]: %s: %f\n", (warn), (f));\
}
#define log_warn_char(warn, c) {\
  printf("[\033[33;1;4mWarning\033[0m]: %s: '%c'\n", (warn), (c));\
}

/**
 *  display error and warning macros, called like this:
 *  log_err(type)(error_meessage, param), eg:
 *  char * filename = "myfile.txt"
 *  log_err(string)("reading file", filename);
 *
 *  types are: none, string, char, int, float
 **/
#define log_err(t) CAT(log_err_, t)
#define log_warn(t) CAT(log_warn_, t)

/**
 *  simple assert macro, display error and exit with exit status 1
 *  if assertion is false
 **/
#define assert(c, msg) if(!(c)) {log_err(none)((msg)); exit(1);}
#define asserts(c, msg, str) if(!(c)) {log_err(string)((msg), (str)); exit(1);}


/**
 *  display warning when assertion is false
 **/
#define assert_warn(c, msg) if(!(c)) {log_warn(none)((msg))}
/**
 *  display warning when assertion is false
 **/
#define assert_err(c, msg) if(!(c)) {log_err(none)((msg))}

#endif /* _ASSERT_H_ */
