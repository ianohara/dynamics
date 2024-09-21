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

static error_t it_takes_two_state(m_t* x, m_t* v, m_t* x_next) {
    m_set_all(x_next, 0);
    m_add(x, v, x_next);

    return E_OK;
}

static error_t direct_measurement_fn(m_t* x, m_t* n, m_t* y) {
    m_set_all(y, 0);
    m_add(x, n, y);

    return E_OK;
}

void test_filtering_kalman__new(void)
{
    m_t* initial_state = m_new(2, 1);
    m_t* initial_state_covariance = m_identity(2);
    m_t* initial_measurement_covariance = m_identity(2);

    kalman_context_t *kt = kalman_new(initial_state, it_takes_two_state, direct_measurement_fn, initial_state_covariance, initial_measurement_covariance);
    cl_assert(kt);
}
