#include <stdint.h>
#include <stdio.h>
#include <math.h>

/* Includes from the testing source tree */
#include "clar.h"
#include "test.h"

/* Includes from the project source tree */
#include "filtering/kalman.h"

void test_filtering_kalman__initialize(void) {
    global_test_counter++;
}

void test_filtering_kalman__cleanup(void)
{
}

// Simple state transition: x_next = x + v (additive input)
static error_t simple_state_fn(m_t* x, m_t* v, m_t* x_next) {
    m_set_all(x_next, 0);
    m_add(x, v, x_next);
    return E_OK;
}

// Measurement function: observe the first state element only
static error_t observe_first_element_fn(m_t* x, m_t* n, m_t* y) {
    (void)n;  // Noise is handled separately in the UKF
    m_set(y, 0, 0, m_get(x, 0, 0));
    return E_OK;
}

void test_filtering_kalman__new(void)
{
    m_t* initial_state = m_new(2, 1);
    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_context_t *kt = kalman_new(initial_state, simple_state_fn, observe_first_element_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    kalman_free(kt);
    m_free(initial_state);
    m_free(process_covariance);
    m_free(measurement_covariance);
}

// Basic smoke test - just verify kalman_step runs without error
void test_filtering_kalman__nonsense_step(void)
{
    m_t* initial_state = m_new(2, 1);
    m_set(initial_state, 0, 0, 0.0);
    m_set(initial_state, 1, 0, 0.0);

    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_context_t *kt = kalman_new(initial_state, simple_state_fn, observe_first_element_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(2, 1);
    m_set(input_vector, 0, 0, 1.0);
    m_set(input_vector, 1, 0, 1.0);

    m_t* measurement = m_new(1, 1);
    m_set(measurement, 0, 0, 0.5);

    error_t result = kalman_step(kt, input_vector, measurement);
    cl_assert_(result == E_OK, "kalman_step should return E_OK");

    kalman_free(kt);
    m_free(initial_state);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
}
