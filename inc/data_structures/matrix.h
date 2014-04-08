#ifndef __MATRIX_H__3434343
#define __MATRIX_H__3434343

#include <stdint.h>

#include "errors.h"

#define MATRIX_DATA_TYPE float

typedef struct m {
    size_t rows;
    size_t cols;
    MATRIX_DATA_TYPE *data;
} m_t;

m_t* m_new(size_t rows, size_t cols);

error_t m_set(m_t *mat, size_t m, size_t n, MATRIX_DATA_TYPE val);

/* Multiple two matricies such that:
     res = lhs*rhs
   It is alright if res is either lhs or rhs as long
   as the sizes work out.
*/
error_t m_mult(m_t *lhs, m_t *rhs, m_t *res);

error_t m_add(m_t *lhs, m_t *rhs, m_t *res);

error_t m_negate(m_t *lhs, m_t *rhs, m_t *res);

#endif /* __MATRIX_H__3434343 */
