#ifndef __VECTOR_H__23232323
#define __VECTOR_H__23232323

#include <stdint.h>

#include "errors.h"

#define VECTOR_DATA_TYPE float

typedef struct v {
    size_t len;
    VECTOR_DATA_TYPE *data;
} v_t;

v_t* v_new(size_t len);

/* Sum the two vs and put the result in res.  It is
   fine for res to be the same as either the lhs or the rhs
*/
error_t v_sum(v_t *lhs, v_t *rhs, v_t *res);

/* Negate a v and put the result in res */
error_t v_negate(v_t *v, v_t *res);

#endif /* __VECTOR_H__23232323 */
