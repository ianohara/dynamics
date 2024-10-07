#include "linear_algebra/decompositions.h"
#include "linear_algebra/properties.h"

error_t la_decompositions_cholesky(m_t* A, m_t* L) {
    if (!A || !L) {
        return E_NULLP;
    }

    if (A->rows != A->cols) {
        return E_VAL;
    }

    if (A->rows != L->rows || A->cols != L->cols) {
        return E_VAL;
    }

    if (!la_is_hermitian(A)) {
        return E_VAL;
    }

    // Algorithm from here: https://en.wikipedia.org/wiki/Cholesky_decomposition#Computation
    /*
    for (i = 0; i < dimensionSize; i++) {
    for (j = 0; j <= i; j++) {
        float sum = 0;
        for (k = 0; k < j; k++)
            sum += L[i][k] * L[j][k];

        if (i == j)
            L[i][j] = sqrt(A[i][i] - sum);
        else
            L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
    }
    }
    */
    if (E_OK != m_set_all(L, 0)) {
        return E_ERR;
    }

    const size_t rows = A->rows;
    for (size_t m = 0; m < rows; m++) {
        for (size_t n = 0; n <= m; n++) {
            m_data_t sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += m_get(L, m, k)*m_get(L, n, k);
            }

            if (m == n) {
                m_set(L, m, n, sqrt(m_get(A, m, m) - sum));
            } else {
                m_set(L, m, n, 1.0 / m_get(L, n, n) * (m_get(A, m, n) - sum));
            }
        }
    }

    return E_OK;
}

// See here for stable gram schmidt: https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
error_t la_decompositions_gram_schmidt(m_t* A, m_t* Q) {
    if (!A || !Q) {
        return E_NULLP;
    }

    if (!m_same_size(A, Q)) {
        return E_VAL;
    }

    if (E_OK != m_set_all(Q, 0)) {
        return E_ERR;
    }

    if (E_OK != m_normalize_column_l2(A, Q, 0)) {
        return E_ERR;
    }

    for (size_t n = 1; n < A->cols; n++) {
        if (E_OK != m_copy_column(A, n, Q, n)) {
            return E_ERR;
        }

        for (size_t n2 = 0; n2 < n; n2++) {
            m_data_t col_dot;
            if (E_OK != m_column_dot_product(Q, n2, Q, n, &col_dot)) {
                return E_ERR;
            }
            for (size_t m = 0; m < Q->rows; m++) {
                m_set(Q, m, n, m_get(Q, m, n) - col_dot*m_get(Q, m, n2));
            }
        }
        if (E_OK != m_normalize_column_l2(Q, Q, n)) {
            return E_ERR;
        }
    }

    return E_OK;
}

/* Uses the Gram-Schmidt process as outlined here: https://en.wikipedia.org/wiki/QR_decomposition#Computing_the_QR_decomposition
 *
 * For now this only supports the square A case.
 */
error_t la_decompositions_qr(m_t* A, m_t* Q, m_t* R) {
    if (!A || !Q || !R) {
        return E_NULLP;
    }

    if (!m_is_square(A) || !m_is_square(Q) || !m_is_square(R)) {
        return E_VAL;
    }

    if (!m_same_size(A, Q) || !m_same_size(A, R)) {
        return E_VAL;
    }

    if (E_OK != m_set_all(R, 0)) {
        return E_ERR;
    }

    // This will zero out Q for us.
    if (E_OK != la_decompositions_gram_schmidt(A, Q)) {
        return E_ERR;
    }

    for (size_t m = 0; m < R->rows; m++) {
        for (size_t n = m; n < R->cols; n++) {
            m_data_t e_m_dot_a_n;
            if (E_OK != m_column_dot_product(Q, m, A, n, &e_m_dot_a_n)) {
                return E_ERR;
            }

            m_set(R, m, n, e_m_dot_a_n);
        }
    }

    return E_OK;
}

error_t la_decompositions_invert_positive_semi_definite(m_t* A, m_t* res) {
    if (!A || !res) {
        return E_NULLP;
    }

    return E_OK;
}