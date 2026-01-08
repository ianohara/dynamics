#include "linear_algebra/properties.h"

bool la_is_hermitian(m_t* A) {
    if (!A || !m_is_square(A)) {
        return false;
    }

    m_t* A_star = m_new(A->cols, A->rows);
    if (!A_star) {
        return false;
    }

    if (E_OK != m_transpose(A, A_star)) {
        m_free(A_star);
        return false;
    }

    // Use epsilon comparison for floating point
    bool result = m_epsilon_equal(A, A_star, 1e-10);
    m_free(A_star);
    return result;
}

