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

/* Invert a positive semi-definite matrix using Cholesky decomposition.
 * A = L * L^T, so A^-1 = (L^-1)^T * L^-1
 */
error_t la_decompositions_invert_positive_semi_definite(m_t* A, m_t* res) {
    if (!A || !res) {
        return E_NULLP;
    }

    if (!m_is_square(A) || !m_same_size(A, res)) {
        return E_VAL;
    }

    size_t n = A->rows;

    // Allocate temporary matrices for L and L_inv
    m_t* L = m_new(n, n);
    m_t* L_inv = m_new(n, n);
    m_t* L_inv_T = m_new(n, n);

    if (!L || !L_inv || !L_inv_T) {
        m_free(L);
        m_free(L_inv);
        m_free(L_inv_T);
        return E_ERR;
    }

    // Compute Cholesky decomposition: A = L * L^T
    if (E_OK != la_decompositions_cholesky(A, L)) {
        m_free(L);
        m_free(L_inv);
        m_free(L_inv_T);
        return E_ERR;
    }

    // Compute L^-1 using forward substitution
    // For each column j of L_inv, solve L * x = e_j
    m_set_all(L_inv, 0);
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) {
            if (i < j) {
                // L_inv[i][j] = 0 (upper triangular part)
                m_set(L_inv, i, j, 0);
            } else if (i == j) {
                // Diagonal: L_inv[i][i] = 1 / L[i][i]
                m_set(L_inv, i, j, 1.0 / m_get(L, i, i));
            } else {
                // Below diagonal: solve using forward substitution
                // L_inv[i][j] = -sum(L[i][k] * L_inv[k][j]) / L[i][i] for k < i
                m_data_t sum = 0;
                for (size_t k = j; k < i; k++) {
                    sum += m_get(L, i, k) * m_get(L_inv, k, j);
                }
                m_set(L_inv, i, j, -sum / m_get(L, i, i));
            }
        }
    }

    // Compute L_inv^T
    m_transpose(L_inv, L_inv_T);

    // Compute A^-1 = L_inv^T * L_inv
    m_mult(L_inv_T, L_inv, res);

    m_free(L);
    m_free(L_inv);
    m_free(L_inv_T);

    return E_OK;
}

/* Cholesky rank-1 update: L*L^T + x*x^T = L'*L'^T
 * Uses the algorithm from LINPACK (dchud).
 * L is modified in place, x is used as workspace and modified.
 */
error_t la_decompositions_cholesky_update(m_t* L, m_t* x) {
    if (!L || !x) {
        return E_NULLP;
    }

    if (!m_is_square(L) || !m_is_vector(x) || L->rows != x->rows) {
        return E_VAL;
    }

    size_t n = L->rows;

    for (size_t k = 0; k < n; k++) {
        m_data_t Lkk = m_get(L, k, k);
        m_data_t xk = m_get(x, k, 0);

        m_data_t r = sqrt(Lkk * Lkk + xk * xk);
        m_data_t c = r / Lkk;
        m_data_t s = xk / Lkk;

        m_set(L, k, k, r);

        for (size_t i = k + 1; i < n; i++) {
            m_data_t Lik = m_get(L, i, k);
            m_data_t xi = m_get(x, i, 0);

            m_set(L, i, k, (Lik + s * xi) / c);
            m_set(x, i, 0, c * xi - s * m_get(L, i, k));
        }
    }

    return E_OK;
}

/* Cholesky rank-1 downdate: L*L^T - x*x^T = L'*L'^T
 * Uses the algorithm from LINPACK (dchdd).
 * L is modified in place, x is used as workspace and modified.
 * Returns E_VAL if result would not be positive definite.
 */
error_t la_decompositions_cholesky_downdate(m_t* L, m_t* x) {
    if (!L || !x) {
        return E_NULLP;
    }

    if (!m_is_square(L) || !m_is_vector(x) || L->rows != x->rows) {
        return E_VAL;
    }

    size_t n = L->rows;

    for (size_t k = 0; k < n; k++) {
        m_data_t Lkk = m_get(L, k, k);
        m_data_t xk = m_get(x, k, 0);

        m_data_t r_sq = Lkk * Lkk - xk * xk;
        if (r_sq <= 0) {
            // Result would not be positive definite
            return E_VAL;
        }

        m_data_t r = sqrt(r_sq);
        m_data_t c = r / Lkk;
        m_data_t s = xk / Lkk;

        m_set(L, k, k, r);

        for (size_t i = k + 1; i < n; i++) {
            m_data_t Lik = m_get(L, i, k);
            m_data_t xi = m_get(x, i, 0);

            m_set(L, i, k, (Lik - s * xi) / c);
            m_set(x, i, 0, c * xi - s * m_get(L, i, k));
        }
    }

    return E_OK;
}