#ifndef __ERRORS_H__121212
#define __ERRORS_H__121212

typedef enum dyn_error {
    E_OK    = 0,
    E_NULLP = 1,
    E_VAL   = 2,
    E_ERR   = 3, // Generic error
    E_INIT  = 4, // Initialization error
} dyn_error_t;

#endif /* __ERRORS_H__121212 */
