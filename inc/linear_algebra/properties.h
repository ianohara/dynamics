#ifndef __LINEAR_ALGEBRA_PROPERTIES__
#define __LINEAR_ALGEBRA_PROPERTIES__

#include <stdbool.h>

#include "data_structures/matrix.h"

bool la_is_hermitian(m_t* A);
bool la_is_positive_definite(m_t* A);

#endif