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

/* Returns the Q and R matricies for the QR decomposition of A.
 *
 * A must be square and real valued.
*/
error_t la_decompositions_qr(m_t* A, m_t* Q, m_t* R);

/* Returns Q where the columns of Q are the orthonormal vectors obtained by carrying out the Gram-Schmidt process on A */
error_t la_decompositions_gram_schmidt(m_t* A, m_t* Q);

error_t la_decompositions_invert_positive_semi_definite(m_t* A, m_t* res);

/* Cholesky rank-1 update: Given L where A = L*L^T, compute L' where A + x*x^T = L'*L'^T
 * L is modified in place. x is a column vector.
 */
error_t la_decompositions_cholesky_update(m_t* L, m_t* x);

/* Cholesky rank-1 downdate: Given L where A = L*L^T, compute L' where A - x*x^T = L'*L'^T
 * L is modified in place. x is a column vector.
 * Returns E_VAL if result would not be positive definite.
 */
error_t la_decompositions_cholesky_downdate(m_t* L, m_t* x);

#endif /* __LINEAR_ALGEBRA_DECOMPOSITONS__ */