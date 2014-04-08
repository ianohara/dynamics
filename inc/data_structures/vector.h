#ifndef __VECTOR_H__23232323
#define __VECTOR_H__23232323

#include <stdint.h>

#include "errors.h"

#define VECTOR_DATA_TYPE float

typedef struct v {
    size_t len;
    VECTOR_DATA_TYPE *data;
} v_t;

/* Returns a new vector of length len.  No values are set */
v_t* v_new(size_t len);

/* Deletes an existing vector and frees any memory it uses */
error_t v_del(v_t *);

/* Sum the two vs and put the result in res.  It is
   fine for res to be the same as either the lhs or the rhs
*/
error_t v_sum(v_t *lhs, v_t *rhs, v_t *res);

/* Negate a v and put the result in res */
error_t v_negate(v_t *v, v_t *res);

#endif /* __VECTOR_H__23232323 */
