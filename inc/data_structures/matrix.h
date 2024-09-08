#ifndef __MATRIX_H__3434343
#define __MATRIX_H__3434343

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "errors.h"

#define M_NAN ((m_data_t)NAN)

typedef double m_data_t;

typedef struct m {
    size_t rows;
    size_t cols;
    m_data_t *data;
} m_t;

m_t* m_new(size_t rows, size_t cols);
error_t m_del(m_t *m);

error_t m_set(m_t *mat, size_t m, size_t n, m_data_t val);
m_data_t m_get(m_t *mat, size_t m, size_t n);

/* Multiple two matricies such that:
     res = lhs*rhs

   It is NOT alright for res to be the lhs or rhs.
*/
error_t m_mult(m_t *lhs, m_t *rhs, m_t *res);

/*  Add two matricies of the same dimensions.

    It is alright for rhs and/or lhs to be the same as res
*/
error_t m_add(m_t *lhs, m_t *rhs, m_t *res);

/*  Negate a matrix.

    It is alright for mat to be the same as res
*/
error_t m_negate(m_t *mat, m_t *res);

/* Get the transpose of the given matrix. */
error_t m_transpose(m_t *mat, m_t *res);

/* Set all entries of the matrix to the given value. */
error_t m_set_all(m_t *mat, m_data_t val);
#endif /* __MATRIX_H__3434343 */
