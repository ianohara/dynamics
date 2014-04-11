#ifndef __VECTOR_H__23232323
#define __VECTOR_H__23232323

#include <stdint.h>
#include <math.h>

#include "errors.h"

typedef float v_data_t;

#define V_NAN ((v_data_t)NAN)

typedef struct v {
    size_t len;
    v_data_t *data;
} v_t;

/* Returns a new vector of length len.  No values are set */
v_t* v_new(size_t len);

/* Returns a new vector with the values in *d filled into it */
v_t* v_new_from_floats(float *d, size_t dlen);

/* Deletes an existing vector and frees any memory it uses */
error_t v_del(v_t *v);

error_t v_set(const v_t *v, size_t ind, v_data_t val);
v_data_t v_get(const v_t * const v, size_t ind);
size_t v_len(const v_t * const v);

/* Sum the two vs and put the result in res.  It is
   fine for res to be the same as either the lhs or the rhs
*/
error_t v_sum(v_t *lhs, v_t *rhs, v_t *res);

/* Negate a v and put the result in res */
error_t v_negate(v_t *v, v_t *res);

#endif /* __VECTOR_H__23232323 */
