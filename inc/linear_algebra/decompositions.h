#ifndef __LINEAR_ALGEBRA_DECOMPOSITONS__
#define __LINEAR_ALGEBRA_DECOMPOSITONS__

#include "errors.h"
#include "data_structures/matrix.h"

/* Returns the lower triangular matrix (L) of the Cholesky decomposition.  The full
 * decomposition is A = LL* where L* is the conjugate transpose of L.
 *
 * A must be Hermitian.  AKA Square and equal to its own conjugate transpose.
 * A must be positive definite.  AKA x^TAx > 0 for all x.
 *
 * See: https://en.wikipedia.org/wiki/Cholesky_decomposition
 */
error_t la_decompositions_cholesky(m_t* A, m_t* L);

#endif /* __LINEAR_ALGEBRA_DECOMPOSITONS__ */