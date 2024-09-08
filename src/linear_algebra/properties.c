#include "linear_algebra/properties.h"

bool la_is_hermitian(m_t* A) {
    m_t* A_star = m_new(A->cols, A->rows);

    if (E_OK != m_transpose(A, A_star)) {
        return false;
    }

    return m_equal(A, A_star);
}

