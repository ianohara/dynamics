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
static dyn_error_t simple_state_fn(m_t* x, m_t* v, m_t* x_next) {
    m_set_all(x_next, 0);
    m_add(x, v, x_next);
    return E_OK;
}

// Measurement function: observe the first state element only
static dyn_error_t observe_first_element_fn(m_t* x, m_t* n, m_t* y) {
    (void)n;  // Noise is handled separately in the UKF
    m_set(y, 0, 0, m_get(x, 0, 0));
    return E_OK;
}

void test_filtering_kalman__new(void)
{
    m_t* initial_state = m_new(2, 1);
    m_t* initial_covariance = m_identity(2);
    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_context_t *kt = kalman_new(initial_state, initial_covariance,
                                       simple_state_fn, observe_first_element_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    kalman_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
}

// Basic smoke test - just verify kalman_step runs without error
void test_filtering_kalman__nonsense_step(void)
{
    m_t* initial_state = m_new(2, 1);
    m_set(initial_state, 0, 0, 0.0);
    m_set(initial_state, 1, 0, 0.0);

    m_t* initial_covariance = m_identity(2);
    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_context_t *kt = kalman_new(initial_state, initial_covariance,
                                       simple_state_fn, observe_first_element_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(2, 1);
    m_set(input_vector, 0, 0, 1.0);
    m_set(input_vector, 1, 0, 1.0);

    m_t* measurement = m_new(1, 1);
    m_set(measurement, 0, 0, 0.5);

    dyn_error_t result = kalman_step(kt, input_vector, measurement);
    cl_assert_(result == E_OK, "kalman_step should return E_OK");

    kalman_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
}

// ============================================================================
// CONSTANT POSITION MODEL TEST
// State: [position]
// The true position is constant, and we observe it with noise.
// The filter should converge to the true position.
// ============================================================================

static dyn_error_t constant_position_state_fn(m_t* x, m_t* v, m_t* x_next) {
    (void)v;  // No input affects constant position
    // Position stays the same
    m_set(x_next, 0, 0, m_get(x, 0, 0));
    return E_OK;
}

static dyn_error_t direct_position_measurement_fn(m_t* x, m_t* n, m_t* y) {
    (void)n;
    // Directly observe position
    m_set(y, 0, 0, m_get(x, 0, 0));
    return E_OK;
}

void test_filtering_kalman__constant_position(void)
{
    // True position is 5.0
    const m_data_t true_position = 5.0;

    // Start with initial guess of 0
    m_t* initial_state = m_new(1, 1);
    m_set(initial_state, 0, 0, 0.0);

    // Large initial uncertainty (we don't know where we are)
    m_t* initial_covariance = m_new(1, 1);
    m_set(initial_covariance, 0, 0, 100.0);

    // Small process noise (position doesn't change)
    m_t* process_covariance = m_new(1, 1);
    m_set(process_covariance, 0, 0, 0.01);

    // Moderate measurement noise
    m_t* measurement_covariance = m_new(1, 1);
    m_set(measurement_covariance, 0, 0, 1.0);

    kalman_context_t *kt = kalman_new(initial_state, initial_covariance,
                                       constant_position_state_fn,
                                       direct_position_measurement_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(1, 1);
    m_set(input_vector, 0, 0, 0.0);

    m_t* measurement = m_new(1, 1);
    m_t* state_out = m_new(1, 1);

    // Simulate noisy measurements around true position
    // Using a simple deterministic "noise" pattern for reproducibility
    m_data_t noise_pattern[] = {0.5, -0.3, 0.8, -0.1, 0.2, -0.6, 0.4, -0.2, 0.1, -0.4};
    size_t num_steps = sizeof(noise_pattern) / sizeof(noise_pattern[0]);

    for (size_t i = 0; i < num_steps; i++) {
        m_set(measurement, 0, 0, true_position + noise_pattern[i]);
        dyn_error_t result = kalman_step(kt, input_vector, measurement);
        cl_assert_(result == E_OK, "kalman_step should succeed");
    }

    // Get final state estimate
    kalman_get_state(kt, state_out);
    m_data_t final_estimate = m_get(state_out, 0, 0);

    // The estimate should be close to the true position (within 1.0)
    m_data_t error = fabs(final_estimate - true_position);
    cl_assert_(error < 1.0, "Final estimate should be close to true position");

    kalman_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
    m_free(state_out);
}

// ============================================================================
// CONSTANT VELOCITY MODEL TEST
// State: [position, velocity]
// The object moves with constant velocity. We only observe position.
// The filter should track both position and velocity.
// ============================================================================

static dyn_error_t constant_velocity_state_fn(m_t* x, m_t* v, m_t* x_next) {
    (void)v;
    // position_next = position + velocity * dt (dt = 1)
    // velocity_next = velocity
    m_data_t position = m_get(x, 0, 0);
    m_data_t velocity = m_get(x, 1, 0);

    m_set(x_next, 0, 0, position + velocity);
    m_set(x_next, 1, 0, velocity);
    return E_OK;
}

static dyn_error_t position_only_measurement_fn(m_t* x, m_t* n, m_t* y) {
    (void)n;
    // Only observe position
    m_set(y, 0, 0, m_get(x, 0, 0));
    return E_OK;
}

void test_filtering_kalman__constant_velocity(void)
{
    // True initial position = 0, velocity = 2
    const m_data_t true_velocity = 2.0;

    // Start with initial guess: position = 0, velocity = 0
    m_t* initial_state = m_new(2, 1);
    m_set(initial_state, 0, 0, 0.0);  // position guess
    m_set(initial_state, 1, 0, 0.0);  // velocity guess (wrong!)

    // Large initial uncertainty
    m_t* initial_covariance = m_identity(2);
    m_scalar_multiply(initial_covariance, 100.0, initial_covariance);

    // Process noise
    m_t* process_covariance = m_identity(2);
    m_scalar_multiply(process_covariance, 0.1, process_covariance);

    // Measurement noise
    m_t* measurement_covariance = m_new(1, 1);
    m_set(measurement_covariance, 0, 0, 1.0);

    kalman_context_t *kt = kalman_new(initial_state, initial_covariance,
                                       constant_velocity_state_fn,
                                       position_only_measurement_fn,
                                       process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(2, 1);
    m_set(input_vector, 0, 0, 0.0);
    m_set(input_vector, 1, 0, 0.0);

    m_t* measurement = m_new(1, 1);
    m_t* state_out = m_new(2, 1);

    // Simulate: true_position(t) = true_velocity * t
    // Add noise pattern for measurements
    m_data_t noise_pattern[] = {0.3, -0.5, 0.2, -0.1, 0.4, -0.3, 0.1, -0.2, 0.5, -0.4,
                                 0.2, -0.3, 0.1, -0.5, 0.3, -0.1, 0.4, -0.2, 0.2, -0.3};
    size_t num_steps = sizeof(noise_pattern) / sizeof(noise_pattern[0]);

    for (size_t t = 0; t < num_steps; t++) {
        m_data_t true_position = true_velocity * (t + 1);
        m_set(measurement, 0, 0, true_position + noise_pattern[t]);

        dyn_error_t result = kalman_step(kt, input_vector, measurement);
        cl_assert_(result == E_OK, "kalman_step should succeed");
    }

    // Get final state estimate
    kalman_get_state(kt, state_out);
    m_data_t final_position = m_get(state_out, 0, 0);
    m_data_t final_velocity = m_get(state_out, 1, 0);

    // Expected final true position
    m_data_t expected_position = true_velocity * num_steps;

    // Position estimate should be reasonable (within 3.0 of true position)
    m_data_t position_error = fabs(final_position - expected_position);
    cl_assert_(position_error < 3.0, "Position estimate should be close to true position");

    // Velocity estimate should have converged toward true velocity (within 1.0)
    m_data_t velocity_error = fabs(final_velocity - true_velocity);
    cl_assert_(velocity_error < 1.0, "Velocity estimate should converge toward true velocity");

    kalman_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
    m_free(state_out);
}

// ============================================================================
// SQUARE-ROOT UKF TESTS
// These tests mirror the standard UKF tests but use the SR-UKF implementation
// ============================================================================

void test_filtering_kalman__sqrt_new(void)
{
    m_t* initial_state = m_new(2, 1);
    m_t* initial_covariance = m_identity(2);
    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_sqrt_context_t *kt = kalman_sqrt_new(initial_state, initial_covariance,
                                                 simple_state_fn, observe_first_element_fn,
                                                 process_covariance, measurement_covariance);
    cl_assert(kt);

    kalman_sqrt_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
}

void test_filtering_kalman__sqrt_nonsense_step(void)
{
    m_t* initial_state = m_new(2, 1);
    m_set(initial_state, 0, 0, 0.0);
    m_set(initial_state, 1, 0, 0.0);

    m_t* initial_covariance = m_identity(2);
    m_t* process_covariance = m_identity(2);
    m_t* measurement_covariance = m_identity(1);

    kalman_sqrt_context_t *kt = kalman_sqrt_new(initial_state, initial_covariance,
                                                 simple_state_fn, observe_first_element_fn,
                                                 process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(2, 1);
    m_set(input_vector, 0, 0, 1.0);
    m_set(input_vector, 1, 0, 1.0);

    m_t* measurement = m_new(1, 1);
    m_set(measurement, 0, 0, 0.5);

    dyn_error_t result = kalman_sqrt_step(kt, input_vector, measurement);
    cl_assert_(result == E_OK, "kalman_sqrt_step should return E_OK");

    kalman_sqrt_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
}

void test_filtering_kalman__sqrt_constant_position(void)
{
    // True position is 5.0
    const m_data_t true_position = 5.0;

    // Start with initial guess of 0
    m_t* initial_state = m_new(1, 1);
    m_set(initial_state, 0, 0, 0.0);

    // Large initial uncertainty (we don't know where we are)
    m_t* initial_covariance = m_new(1, 1);
    m_set(initial_covariance, 0, 0, 100.0);

    // Small process noise (position doesn't change)
    m_t* process_covariance = m_new(1, 1);
    m_set(process_covariance, 0, 0, 0.01);

    // Moderate measurement noise
    m_t* measurement_covariance = m_new(1, 1);
    m_set(measurement_covariance, 0, 0, 1.0);

    kalman_sqrt_context_t *kt = kalman_sqrt_new(initial_state, initial_covariance,
                                                 constant_position_state_fn,
                                                 direct_position_measurement_fn,
                                                 process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(1, 1);
    m_set(input_vector, 0, 0, 0.0);

    m_t* measurement = m_new(1, 1);
    m_t* state_out = m_new(1, 1);

    // Simulate noisy measurements around true position
    m_data_t noise_pattern[] = {0.5, -0.3, 0.8, -0.1, 0.2, -0.6, 0.4, -0.2, 0.1, -0.4};
    size_t num_steps = sizeof(noise_pattern) / sizeof(noise_pattern[0]);

    for (size_t i = 0; i < num_steps; i++) {
        m_set(measurement, 0, 0, true_position + noise_pattern[i]);
        dyn_error_t result = kalman_sqrt_step(kt, input_vector, measurement);
        cl_assert_(result == E_OK, "kalman_sqrt_step should succeed");
    }

    // Get final state estimate
    kalman_sqrt_get_state(kt, state_out);
    m_data_t final_estimate = m_get(state_out, 0, 0);

    // The estimate should be close to the true position (within 1.0)
    m_data_t error = fabs(final_estimate - true_position);
    cl_assert_(error < 1.0, "Final estimate should be close to true position");

    kalman_sqrt_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
    m_free(state_out);
}

void test_filtering_kalman__sqrt_constant_velocity(void)
{
    // True initial position = 0, velocity = 2
    const m_data_t true_velocity = 2.0;

    // Start with initial guess: position = 0, velocity = 0
    m_t* initial_state = m_new(2, 1);
    m_set(initial_state, 0, 0, 0.0);  // position guess
    m_set(initial_state, 1, 0, 0.0);  // velocity guess (wrong!)

    // Large initial uncertainty
    m_t* initial_covariance = m_identity(2);
    m_scalar_multiply(initial_covariance, 100.0, initial_covariance);

    // Process noise
    m_t* process_covariance = m_identity(2);
    m_scalar_multiply(process_covariance, 0.1, process_covariance);

    // Measurement noise
    m_t* measurement_covariance = m_new(1, 1);
    m_set(measurement_covariance, 0, 0, 1.0);

    kalman_sqrt_context_t *kt = kalman_sqrt_new(initial_state, initial_covariance,
                                                 constant_velocity_state_fn,
                                                 position_only_measurement_fn,
                                                 process_covariance, measurement_covariance);
    cl_assert(kt);

    m_t* input_vector = m_new(2, 1);
    m_set(input_vector, 0, 0, 0.0);
    m_set(input_vector, 1, 0, 0.0);

    m_t* measurement = m_new(1, 1);
    m_t* state_out = m_new(2, 1);

    // Simulate: true_position(t) = true_velocity * t
    m_data_t noise_pattern[] = {0.3, -0.5, 0.2, -0.1, 0.4, -0.3, 0.1, -0.2, 0.5, -0.4,
                                 0.2, -0.3, 0.1, -0.5, 0.3, -0.1, 0.4, -0.2, 0.2, -0.3};
    size_t num_steps = sizeof(noise_pattern) / sizeof(noise_pattern[0]);

    for (size_t t = 0; t < num_steps; t++) {
        m_data_t true_position = true_velocity * (t + 1);
        m_set(measurement, 0, 0, true_position + noise_pattern[t]);

        dyn_error_t result = kalman_sqrt_step(kt, input_vector, measurement);
        cl_assert_(result == E_OK, "kalman_sqrt_step should succeed");
    }

    // Get final state estimate
    kalman_sqrt_get_state(kt, state_out);
    m_data_t final_position = m_get(state_out, 0, 0);
    m_data_t final_velocity = m_get(state_out, 1, 0);

    // Expected final true position
    m_data_t expected_position = true_velocity * num_steps;

    // Position estimate should be reasonable (within 3.0 of true position)
    m_data_t position_error = fabs(final_position - expected_position);
    cl_assert_(position_error < 3.0, "Position estimate should be close to true position");

    // Velocity estimate should have converged toward true velocity (within 1.0)
    m_data_t velocity_error = fabs(final_velocity - true_velocity);
    cl_assert_(velocity_error < 1.0, "Velocity estimate should converge toward true velocity");

    kalman_sqrt_free(kt);
    m_free(initial_state);
    m_free(initial_covariance);
    m_free(process_covariance);
    m_free(measurement_covariance);
    m_free(input_vector);
    m_free(measurement);
    m_free(state_out);
}
