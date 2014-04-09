#include <stdlib.h>

#include "data_structures/matrix.h"

m_t* m_new(size_t rows, size_t cols)
{
    m_t *nm = NULL;
    m_data_t *d = NULL;
    if (rows <= 0) goto fail;
    if (cols <= 0) goto fail;

    nm = malloc(sizeof *nm);
    if (!nm) goto fail_nm;

    d = malloc(rows*cols*sizeof *d);
    if (!d) goto fail_d;

    nm->rows = rows;
    nm->cols = cols;
    nm->data = d;
    goto out;

    fail_d:
    free(nm);

    fail_nm:
    fail:
    nm = NULL;

    out:
    return nm;
}

error_t m_del(m_t *m)
{
    if (!m) return E_OK;
    if (m->data) free(m->data);
    free(m);
    return E_OK;
}

/* Unsafe - does no checks.  Make sure you have your stuff right!

    For now store as rows because this makes sense for multiplying
    by vectors.
*/
static inline size_t m_get_index(m_t *mat, size_t m, size_t n)
{
    return mat->data[m*mat->cols+n];
}

error_t m_set(m_t *mat, size_t m, size_t n, m_data_t val)
{
    if (!mat) return E_NULLP;
    if (m >= mat->rows) return E_VAL;
    if (n >= mat->cols) return E_VAL;

    mat->data[m_get_index(mat,m,n)] = val;
    return E_OK;
}

m_data_t m_get(m_t *mat, size_t m, size_t n)
{
    if (!mat) return M_NAN;
    if (m >= mat->rows) return M_NAN;
    if (n >= mat->cols) return M_NAN;

    return mat->data[m_get_index(mat,m,n)];
}

error_t m_mult(m_t *lhs, m_t *rhs, m_t *res)
{
    size_t common_dim;
    m_data_t rc_sum;

    /* TODO: Better errors.  4 different failures give E_VAL! */
    if (!lhs || !rhs || !res) return E_NULLP;
    if (lhs->cols != rhs->rows) return E_VAL;
    if (lhs->rows != res->rows) return E_VAL;
    if (rhs->cols != res->cols) return E_VAL;
    if (res == lhs || res == rhs) return E_VAL;

    common_dim = lhs->cols;

    for (size_t m=0; m < lhs->rows; m++) {
        for (size_t n=0; n < rhs->cols; n++) {
            rc_sum = 0.0;
            for (size_t i=0; i < common_dim; i++) {
                rc_sum += m_get(lhs, m, i)*m_get(rhs, i, n);
            }
            m_set(res, m, n, rc_sum);
        }
    }
    return E_OK;
}

error_t m_add(m_t *lhs, m_t *rhs, m_t *res)
{
    if (!lhs || !rhs || !res) return E_NULLP;
    if (lhs->cols != rhs->cols) return E_VAL;
    if (lhs->rows != rhs->rows) return E_VAL;
    if (lhs->cols != res->cols) return E_VAL;
    if (lhs->rows != res->rows) return E_VAL;

    for (size_t m=0; m < res->rows; m++) {
        for (size_t n=0; n < res->cols; n++) {
            m_set(res, m, n, m_get(lhs,m,n)+m_get(rhs,m,n));
        }
    }
    return E_OK;
}

