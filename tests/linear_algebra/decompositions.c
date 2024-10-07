#include <stdint.h>
#include <stdio.h>
#include <math.h>

/* Includes from the testing source tree */
#include "clar.h"
#include "test.h"

/* Includes from the project source tree */
#include "data_structures/matrix.h"
#include "linear_algebra/decompositions.h"

void test_linear_algebra_decompositions__initialize(void)
{
    global_test_counter++;
}

void test_linear_algebra_decompositions__cleanup(void)
{
}

void test_linear_algebra_decompositions__graham_schmidt_basic_cases(void)
{
    // Identity maps to identity
    {
        m_t* identity = m_identity(6);
        m_t* result = m_new(6, 6);
        m_set_all(result, 0);

        cl_assert(la_decompositions_gram_schmidt(identity, result) == E_OK);
        cl_assert(m_epsilon_equal(identity, result, ELEMENT_EPSILON));

        m_free(identity);
        m_free(result);
    }
    // Scaled identity maps to identity
    {
        m_t* identity = m_identity(6);
        m_t* scaling = m_new(6, 6);
        for (size_t n = 0; n < scaling->cols; n++) {
            m_set(scaling, n, n, n + 2);
        }
        m_t* scaled_identity = m_new(6, 6);
        cl_assert(m_mult(identity, scaling, scaled_identity) == E_OK);

        m_t* result = m_new(6, 6);
        m_set_all(result, 0);

        cl_assert(la_decompositions_gram_schmidt(scaled_identity, result) == E_OK);
        cl_assert(m_epsilon_equal(identity, result, ELEMENT_EPSILON));

        m_free(identity);
        m_free(result);
        m_free(scaling);
        m_free(scaled_identity);
    }
    // Already orthogonal maps to same
    {
        double orthogonal_data[36] = {
            -0.39440783, -0.43585156,  0.13539159,  0.7405485 , -0.29391025,
            0.03673812, -0.35429019, -0.54417173, -0.01128728, -0.37183873,
            0.40588416,  0.52461624, -0.4964513 ,  0.44413828, -0.50517339,
            0.27844161,  0.46997312, -0.05169805, -0.46934734,  0.19918369,
            0.73963638, -0.17954912,  0.17786052, -0.35931145, -0.31428041,
            -0.33963319, -0.42209151, -0.38127933, -0.30564933, -0.60738821,
            -0.39166049,  0.40202847, -0.03344176, -0.24121014, -0.63477886,
            0.47194309
        };

        m_t* orthogonal = m_new_from(6, 6, orthogonal_data);
        m_t* results = m_new(6, 6);
        m_set_all(results, 0);

        cl_assert(la_decompositions_gram_schmidt(orthogonal, results) == E_OK);
        // Only use 1e-8 for the epsilon since the input data is only to about 1e-8
        // (just whatever ipython natively printed out).
        cl_assert(m_epsilon_equal(orthogonal, results, 1e-8));

        m_free(orthogonal);
        m_free(results);
    }

    // wikipedia example (from the QR decomposition page): https://en.wikipedia.org/wiki/QR_decomposition
    {
        double A_data[9] = {12.0, -51.0, 4.0, 6.0, 167.0, -68.0, -4.0, 24.0, -41.0};
        double expected_data[9] = {6/7.0, -69/175.0, -58/175.0, 3.0/7.0, 158.0/175.0, 6/175.0, -2.0/7.0, 6/35.0, -33.0/35.0};

        m_t* A = m_new_from(3, 3, A_data);
        m_t* expected = m_new_from(3, 3, expected_data);

        m_t* result = m_new(3, 3);
        m_set_all(result, 0);

        cl_assert(la_decompositions_gram_schmidt(A, result) == E_OK);
        cl_assert(m_epsilon_equal(expected, result, ELEMENT_EPSILON));

        m_free(A);
        m_free(expected);
        m_free(result);
    }

    // smooth failure in invalid input case
}