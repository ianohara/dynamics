#ifndef __MATRIX_H__3434343
#define __MATRIX_H__3434343

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "errors.h"

#define M_NAN ((m_data_t)NAN)

typedef double m_data_t;

typedef struct m
{
    size_t rows;
    size_t cols;
    m_data_t *data;
} m_t;

m_t *m_new(size_t rows, size_t cols);
error_t m_free(m_t *m);

error_t m_set(m_t *mat, size_t m, size_t n, m_data_t val);
m_data_t m_get(m_t *mat, size_t m, size_t n);

/* Multiple two matricies such that:
     res = lhs*rhs

   It is NOT alright for res to be the lhs or rhs.
*/
error_t m_mult(m_t *lhs, m_t *rhs, m_t *res);

/*  Add two matricies of the same dimensions.

    It is alright for rhs and/or lhs to be the same as res
*/
error_t m_add(m_t *lhs, m_t *rhs, m_t *res);

/*  Negate a matrix.

    It is alright for mat to be the same as res
*/
error_t m_negate(m_t *mat, m_t *res);

/* Get the transpose of the given matrix. */
error_t m_transpose(m_t *mat, m_t *res);

/* Returns true if both matricies have exactly equal sizes and components */
bool m_equal(m_t *a, m_t *b);

/* Set all entries of the matrix to the given value. */
error_t m_set_all(m_t *mat, m_data_t val);

/* Returns true if mat is a square matrix. */
bool m_is_square(m_t *mat);

/* Returns true if mats a and b have exactly equal rows and column counts. */
bool m_same_size(m_t *a, m_t *b);

/* Calculates the normalized vector of column col_idx (using the l2 norm) in src and places it in dest.  It is okay if dest == src. */
error_t m_normalize_column_l2(m_t *src, m_t *dest, size_t col_idx);

/* Copy column col from src to dest */
error_t m_copy_column(m_t* src, size_t src_col, m_t* dest, size_t dest_col);

/* Copy the given matrix into the specified position in the destination matrix. */
error_t m_copy_into(m_t* src, m_t* dest, size_t dest_row, size_t dest_col);

/* Add the column from src, scaled by this factor, to the column in dest. */
error_t m_add_scaled_column(m_t* src, size_t src_col_idx, m_data_t src_scale, m_t* dest, size_t dest_col_idx);

/* Calculate the dot product of the two requested columns in A and B and put the result in res. */
error_t m_column_dot_product(m_t* A, size_t a_col, m_t* B, size_t b_col, m_data_t* res);

bool m_is_vector(m_t* mat);

size_t m_max_dim(m_t* mat);

error_t m_outer_product(m_t* lhs, m_t* rhs, m_t* res);
#endif /* __MATRIX_H__3434343 */
