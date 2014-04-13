#ifndef __VECTOR_H__23232323
#define __VECTOR_H__23232323

#include <stdint.h>
#include <math.h>

#include "errors.h"

typedef double v_data_t;

#define V_NAN ((v_data_t)NAN)
#define v_isnan(val) (isnan(val))

/* Iterate over a vector from first index to last
   using the index variable "ind". You are responsible
   for brackets if you want 'em!
*/
#define v_foreach(v,ind) \
    for (size_t __vlen=v_len(v), ind=0; \
         ind < __vlen; \
         ind++)

typedef struct v {
    size_t len;
    v_data_t *data;
} v_t;

/* Returns a new vector of length len.  No values are set */
v_t* v_new(size_t len);
/* Deletes an existing vector and frees any memory it uses */
error_t v_del(v_t *v);

/****
 * A number of helpful vector creating functions.  These all
 * return a vector to a new vector that later needs to be
 * freed with v_del by the caller.
 *****/

/* Returns a new vector with the values in *d filled into it */
v_t* v_new_from_floats(float *d, size_t dlen);
v_t* v_new_from_value(v_data_t val, size_t len);
v_t* v_new_zeros(size_t len);
v_t* v_new_ones(size_t len);

/****
 * Functions for modifying and getting properties of
 * vectors.
 ****/
error_t v_set(const v_t *v, size_t ind, v_data_t val);
v_data_t v_get(const v_t * const v, size_t ind);
size_t v_len(const v_t * const v);


/****
 * Functions for carrying out calculations on
 * vectors.  For all of these, you are required
 * to pass in an existing vector in which the result
 * should be stored.  Only in specifically noted cases
 * can this result vector pointer be the same as one
 * of the operands.
 ****/

/* Multiplication by a scalar.

    It is fine for res to be the same as v.
 */
error_t v_sp(v_data_t s, v_t *v, v_t *res);

/* Sum the two vs and put the result in res.

   It is fine for res to be the same as either the lhs or the rhs
*/
error_t v_sum(v_t *lhs, v_t *rhs, v_t *res);

/* Negate a v and put the result in res

    It is fine for res to be the same as v.
*/
error_t v_negate(v_t *v, v_t *res);

/* Return the dot product of lhs and rhs.
 *
 * Returns V_NAN if there was an error or if any
 * of the values in the vectors to be dotted are V_NAN.
 */
v_data_t v_dot(v_t *lhs, v_t *rhs);

#endif /* __VECTOR_H__23232323 */
