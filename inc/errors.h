#ifndef __ERRORS_H__121212
#define __ERRORS_H__121212

#ifdef error
#   error "error already defined! Use something other than error for the error enum"
#endif /* error */

typedef enum error {
    E_OK    = 0,
    E_NULLP = 1,
    E_VAL   = 2,
} error_t;

#endif /* __ERRORS_H__121212 */
