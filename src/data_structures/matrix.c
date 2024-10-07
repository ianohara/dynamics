#include <stdlib.h>

#include "data_structures/matrix.h"

m_t* m_new(size_t rows, size_t cols)
{
    m_t *nm = NULL;
    m_data_t *d = NULL;

    if (rows <= 0) {
        m_free(nm);
        return NULL;
    }

    if (cols <= 0) {
        m_free(nm);
        return NULL;
    }

    nm = malloc(sizeof *nm);
    if (!nm) {
        m_free(nm);
        return NULL;
    }

    d = malloc(rows*cols*sizeof *d);
    if (!d) {
        m_free(nm);
        return NULL;
    }

    nm->rows = rows;
    nm->cols = cols;
    nm->data = d;

    return nm;
}

m_t* m_new_from(size_t rows, size_t cols, m_data_t* data) {
    m_t* res = m_new(rows, cols);

    if (!res) {
        return NULL;
    }

    size_t data_idx = 0;
    for (size_t m = 0; m < rows; m++) {
        for (size_t n = 0; n < cols; n++) {
            m_set(res, m, n, data[data_idx++]);
        }
    }

    return res;
}

m_t* m_identity(size_t dim) {
    m_t* m = m_new(dim, dim);

    if (!m) {
        return NULL;
    }

    m_set_all(m, 0);
    for (size_t d = 0; d < dim; d++) {
        m_set(m, d, d, 1);
    }

    return m;
}

error_t m_free(m_t *m)
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
    return m*mat->cols+n;
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

error_t m_negate(m_t *mat, m_t *res)
{
    if (!mat || !res) return E_NULLP;
    if (mat->cols != res->cols) return E_VAL;
    if (mat->rows != res->rows) return E_VAL;

    for (size_t m=0; m < res->rows; m++) {
        for (size_t n=0; n < res->cols; n++) {
            m_set(res, m, n, -m_get(mat,m,n));
        }
    }

    return E_OK;
}

error_t m_transpose(m_t *mat, m_t *res) {
    if (!mat || !res) return E_NULLP;

    if (mat->rows != res->cols) {
        return E_VAL;
    }
    if (mat->cols != res->rows) {
        return E_VAL;
    }

    for (size_t m = 0; m < mat->rows; m++) {
        for (size_t n = 0; n < mat->cols; n++) {
            m_set(res, n, m, m_get(mat, m, n));
        }
    }

    return E_OK;
}

bool m_equal(m_t *a, m_t *b) {
    if (!a || !b) {
        return false;
    }

    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }

    for (size_t m = 0; m < a->rows; m++) {
        for (size_t n = 0; n < a->cols; n++) {
            if (m_get(a, m, n) != m_get(b, m, n)) {
                return false;
            }
        }
    }

    return true;
}

bool m_epsilon_equal(m_t *a, m_t *b, m_data_t epsilon) {
    if (!a || !b) {
        return false;
    }

    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }

    for (size_t m = 0; m < a->rows; m++) {
        for (size_t n = 0; n < a->cols; n++) {
            if (fabs(m_get(a, m, n) - m_get(b, m, n)) >= epsilon) {
                return false;
            }
        }
    }

    return true;
}

error_t m_set_all(m_t* mat, m_data_t val) {
    if (!mat) {
        return E_NULLP;
    }

    for (size_t m = 0; m < mat->rows; m++) {
        for (size_t n = 0; n < mat->cols; n++) {
            m_set(mat, m, n, val);
        }
    }

    return E_OK;
}

bool m_is_square(m_t* mat) {
    return mat && mat->rows == mat->cols;
}

bool m_same_size(m_t* a, m_t* b) {
    return a && b && a->rows == b->rows && a->cols == b->cols;
}

error_t m_normalize_column_l2(m_t* src, m_t* dest, size_t col_idx) {
    if (!src || !dest) {
        return E_NULLP;
    }

    if (src->rows != dest->cols || src->cols - 1 < col_idx || dest->cols - 1 < col_idx) {
        return E_VAL;
    }

    m_data_t col_length = 0.0;
    for (size_t m = 0; m < src->rows; m++) {
        col_length += m_get(src, m, col_idx) * m_get(src, m, col_idx);
    }
    col_length = sqrt(col_length);

    for (size_t m = 0; m < src->rows; m++) {
        m_set(dest, m, col_idx, m_get(src, m, col_idx) / col_length);
    }

    return E_OK;
}

error_t m_copy_column(m_t* src, size_t src_col, m_t* dest, size_t dest_col) {
    if (!src || !dest){
        return E_NULLP;
    }

    if (src->rows != dest->rows) {
        return E_VAL;
    }

    if (src_col > src->cols - 1 || dest_col > dest->cols - 1) {
        return E_VAL;
    }

    for (size_t m = 0; m < src->rows; m++) {
        m_set(dest, m, dest_col, m_get(src, m, src_col));
    }

    return E_OK;
}

error_t m_copy_into(m_t* src, m_t* dest, size_t dest_row, size_t dest_col) {
    if (!src || !dest) {
        return E_NULLP;
    }

    if (dest_row + src->rows > dest->rows - 1) {
        return E_VAL;
    }

    if (dest_col + src->cols > dest->cols - 1) {
        return E_VAL;
    }

    for (size_t m = 0; m < src->rows; m++) {
        for (size_t n = 0; n < src->cols; n++) {
            m_set(dest, dest_row + m, dest_col + n, m_get(src, m, n));
        }
    }

    return E_OK;
}

error_t m_add_scaled_column(m_t* src, size_t src_col, m_data_t scale, m_t* dest, size_t dest_col) {
    if (!src || !dest){
        return E_NULLP;
    }

    if (src->rows != dest->rows) {
        return E_VAL;
    }

    if (src_col > src->cols - 1 || dest_col > dest->cols - 1) {
        return E_VAL;
    }

    for (size_t m = 0; m < src->rows; m++) {
        m_set(dest, m, dest_col, m_get(dest, m, dest_col) + scale*m_get(src, m, src_col));
    }

    return E_OK;
}

error_t m_column_dot_product(m_t* A, size_t a_col, m_t* B, size_t b_col, m_data_t* res) {
    if (!A || !B || !res) {
        return E_NULLP;
    }

    if (!m_same_size(A, B)) {
        return E_VAL;
    }

    if (a_col > A->cols - 1 || b_col > B->cols) {
        return E_VAL;
    }

    m_data_t tmp_res = 0.0;

    for (size_t m = 0; m < A->rows; m++) {
        tmp_res += m_get(A, m, a_col)*m_get(B, m, b_col);
    }

    *res = tmp_res;

    return E_OK;
}

bool m_is_vector(m_t* mat) {
    return mat && mat->cols == 1 && mat->rows > 0;
}

size_t m_max_dim(m_t* mat) {
    if (!mat) {
        return 0;
    }

    return mat->rows > mat->cols ? mat->rows : mat->cols;
}

error_t m_outer_product(m_t* lhs, m_t* rhs, m_t* res) {
    if (!lhs || !rhs || !res) {
        return E_NULLP;
    }

    if (!m_is_vector(lhs) || !m_is_vector(rhs)) {
        return E_VAL;
    }

    if (res->rows != lhs->rows || res->cols != rhs->rows) {
        return E_VAL;
    }

    for (size_t m_l = 0; m_l < lhs->rows; m_l++) {
        for (size_t m_r = 0; m_r < rhs->rows; m_r++) {
            m_set(res, m_l, m_r, m_get(lhs, m_l, 0)*m_get(rhs, m_r, 0));
        }
    }

    return E_OK;
}

error_t m_scalar_multiply(m_t* src, m_data_t multiplier, m_t* dest) {
    if (!src || !dest) {
        return E_NULLP;
    }

    if (!m_same_size(src, dest)) {
        return E_VAL;
    }

    for (size_t m = 0; m < src->rows; m++) {
        for (size_t n = 0; n < src->cols; n++) {
            m_set(dest, m, n, multiplier * m_get(src, m, n));
        }
    }

    return E_OK;
}