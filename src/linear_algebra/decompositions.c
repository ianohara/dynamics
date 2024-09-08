#include "linear_algebra/decompositions.h"

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