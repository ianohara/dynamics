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

inline error_t v_set(const v_t *v, size_t ind, v_data_t val) {
    if (!v) return E_NULLP;
    if (ind >= v->len) return E_VAL;
    v->data[ind] = val;
    return E_OK;
}

inline v_data_t v_get(const v_t * const v, size_t ind) {
    if (!v) return V_NAN;
    if (ind >= v->len) return V_NAN;
    return v->data[ind];
}

inline size_t v_len(const v_t * const v) {
    if (!v) return 0;
    return v->len;
}

/* Returns a new vector of length len.  No values are set */
v_t* v_new(size_t len);

/* Deletes an existing vector and frees any memory it uses */
error_t v_del(v_t *v);

/* Sum the two vs and put the result in res.  It is
   fine for res to be the same as either the lhs or the rhs
*/
error_t v_sum(v_t *lhs, v_t *rhs, v_t *res);

/* Negate a v and put the result in res */
error_t v_negate(v_t *v, v_t *res);

#endif /* __VECTOR_H__23232323 */
