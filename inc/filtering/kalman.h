#ifndef __FILTERING_KALMAN_H__
#define __FILTERING_KALMAN_H__

#include "errors.h"
#include "data_structures/matrix.h"

typedef error_t (*kalman_state_fn)(m_t *x, m_t *v, m_t *x_next);
typedef error_t (*kalman_state_to_measurement_fn)(m_t* x, m_t *n, m_t *y_next);
// This uses the notation from "The Unscented Kalman Filter for Nonlinear Estimation" by Wan and Merwe
// The implementation of kalman_step is from the same.
typedef struct kalman_context
{
    /* Alpha must hold to: 1e-4 <= alpha <= 1*/
    m_data_t alpha;

    /* Beta of 2 is best for gaussian distributed state vectors */
    m_data_t beta;

    /* Kappa is normally set to 0 */
    m_data_t kappa;

    // Dimensions
    size_t state_len;
    size_t measurement_len;

    // Any of the fields with leading underscores are internal scratch pad values that you
    // should not touch.
    m_data_t lambda;
    m_data_t eta; // Not mentioned in the paper, but sqrt(state_len + lambda) is used a ton.

    m_t *sigma_weights_m;
    m_t *sigma_weights_c;

    // Process and measurement noise covariances
    m_t *Q;  // Process noise covariance
    m_t *R;  // Measurement noise covariance

    m_t *P_km1;
    m_t *sqrt_P_km1;

    m_t *P_k_t;
    m_t *P_yy;
    m_t *P_xy;

    m_t *x_hat;
    m_t *x_hatm1;
    m_t* x_hat_k_t;
    m_t *chi_k;
    m_t *chi_km1;

    // Various internal scratch matricies
    // Used as the output for the state function
    m_t* x_b;
    m_t* x_a;
    m_t *P_a;
    m_t *P_b;  // Additional scratch for covariance accumulation

    m_t *Y_k;
    m_t *y_hat_k_t;
    m_t *y_a;
    m_t *y_b;  // Additional scratch for measurement covariance accumulation

    m_t *K;
    m_t *K_transpose;
    m_t *scratch_yx;

    kalman_state_fn state_fn;
    kalman_state_to_measurement_fn measurement_fn;
} kalman_context_t;

/* Create a new kalman context.  This does not take ownership (or use - it's fine to immediately free them) of any of the m_t* passed in, so you must manage those after the call.
 *
 * initial_covariance: Initial state covariance estimate. If NULL, uses process_covariance as default.
 */
kalman_context_t *kalman_new(
    m_t *initial_state_guess,
    m_t *initial_covariance,
    kalman_state_fn state_fn,
    kalman_state_to_measurement_fn measurement_fn,
    m_t* process_covariance,
    m_t* measurement_covariance);
void kalman_free(kalman_context_t *context);

error_t kalman_step(kalman_context_t* context, m_t* input_vector, m_t* measurement);

/* Get the current state estimate. Copies the internal state into the provided matrix. */
error_t kalman_get_state(kalman_context_t* context, m_t* state_out);

/* Get the current state covariance estimate. Copies the internal covariance into the provided matrix. */
error_t kalman_get_covariance(kalman_context_t* context, m_t* covariance_out);

#endif