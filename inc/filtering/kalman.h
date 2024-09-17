#ifndef __FILTERING_KALMAN_H__
#define __FILTERING_KALMAN_H__

#include "errors.h"
#include "data_structures/matrix.h"

typedef error_t (*kalman_state_fn)(m_t *x, m_t *u, m_t *x_next);

// This uses the notation from "The square-root unscented Kalman filter..." by Van Der Merwe and Wan, 2001.
typedef struct kalman_context
{
    /* Alpha must hold to: 1e-4 <= alpha <= 1*/
    double alpha;

    /* Beta of 2 is best for gaussian distributed state vectors */
    double beta;

    // Any of the fields with leading underscores are internal scratch pad values that you
    // should not touch.
    double lambda;
    double eta;
    m_t *sigma_weights;

    m_t *x_hat;
    m_t *x_hatm1;
    m_t* x_hat_k_t;
    m_t *S_k;
    m_t *S_k_t;
    m_t *chi_k;
    m_t *chi_km1;

    // Various internal scratch matricies
    // Used as the output for the state function
    m_t* x_scratch;
    m_t* weighted_chi_minus_x_k_t;

    m_t* sqrt_state_covariance;
    m_t* sqrt_measurement_covariance;

    m_t* state_qr_block;
    m_t* state_qr_Q_placeholder;

    m_t* measurement_qr_block;
    m_t* measurement_qr_Q_placeholder;

    kalman_state_fn state_fn;
} kalman_context_t;

/* Create a new kalman context.  This does not take ownership (or use - it's fine to immediately free them) of any of the m_t* passed in, so you must manage those after the call. */
kalman_context_t *kalman_new(m_t *initial_state_guess, kalman_state_fn state_fn, m_t* state_covariance, m_t* measurement_covariance);
void kalman_free(kalman_context_t *context);
#endif