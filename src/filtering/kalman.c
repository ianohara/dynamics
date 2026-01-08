#include <string.h>

#include "filtering/kalman.h"
#include "linear_algebra/decompositions.h"
#include "data_structures/matrix.h"

kalman_context_t* kalman_new(
    m_t *initial_state_guess,
    kalman_state_fn state_fn,
    kalman_state_to_measurement_fn measurement_fn,
    m_t* process_covariance,
    m_t* measurement_covariance) {
    if (!initial_state_guess || !state_fn || !process_covariance || !measurement_covariance) {
        return NULL;
    }

    if (!m_is_vector(initial_state_guess)) {
        return NULL;
    }

    if (!m_is_square(process_covariance) || process_covariance->rows != initial_state_guess->rows) {
        return NULL;
    }

    if (!m_is_square(measurement_covariance)) {
        return NULL;
    }

    size_t state_len = initial_state_guess->rows;
    size_t measurement_len = measurement_covariance->rows;

    kalman_context_t* context = malloc(sizeof(*context));
    if (!context) {
        return NULL;
    }

    // Zero out context to make cleanup safer
    memset(context, 0, sizeof(*context));

    context->state_len = state_len;
    context->measurement_len = measurement_len;

    context->alpha = 1e-3;
    context->beta = 2.0;
    context->kappa = 0.0;

    context->lambda = context->alpha * context->alpha * (state_len + context->kappa) - state_len;
    context->eta = sqrt(state_len + context->lambda);

    context->state_fn = state_fn;
    context->measurement_fn = measurement_fn;

    size_t sigma_point_count = 2 * state_len + 1;

    // Allocate sigma weights
    if (!(context->sigma_weights_m = m_new(sigma_point_count, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->sigma_weights_c = m_new(sigma_point_count, 1))) {
        kalman_free(context);
        return NULL;
    }

    // Calculate weights
    m_data_t weight_0_m = context->lambda / (state_len + context->lambda);
    m_data_t weight_0_c = weight_0_m + (1 - context->alpha * context->alpha + context->beta);
    m_data_t weight_i = 1.0 / (2.0 * (state_len + context->lambda));

    m_set(context->sigma_weights_m, 0, 0, weight_0_m);
    m_set(context->sigma_weights_c, 0, 0, weight_0_c);
    for (size_t i = 1; i < sigma_point_count; i++) {
        m_set(context->sigma_weights_m, i, 0, weight_i);
        m_set(context->sigma_weights_c, i, 0, weight_i);
    }

    // Allocate Q and R (copies of the input covariances)
    if (!(context->Q = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }
    m_copy(process_covariance, context->Q);

    if (!(context->R = m_new(measurement_len, measurement_len))) {
        kalman_free(context);
        return NULL;
    }
    m_copy(measurement_covariance, context->R);

    // Allocate covariance matrices
    if (!(context->P_km1 = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->sqrt_P_km1 = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->P_k_t = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->P_yy = m_new(measurement_len, measurement_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->P_xy = m_new(state_len, measurement_len))) {
        kalman_free(context);
        return NULL;
    }

    // Allocate state vectors
    if (!(context->x_hat = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->x_hatm1 = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->x_hat_k_t = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    // Allocate sigma point matrices
    if (!(context->chi_k = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->chi_km1 = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    // Allocate scratch matrices for state
    if (!(context->x_a = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->x_b = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->P_a = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->P_b = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    // Allocate measurement matrices
    if (!(context->Y_k = m_new(measurement_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->y_hat_k_t = m_new(measurement_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->y_a = m_new(measurement_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->y_b = m_new(measurement_len, measurement_len))) {
        kalman_free(context);
        return NULL;
    }

    // Allocate Kalman gain matrices
    if (!(context->K = m_new(state_len, measurement_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->K_transpose = m_new(measurement_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->scratch_yx = m_new(measurement_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    // Initialize state estimate from guess
    m_copy(initial_state_guess, context->x_hat);

    // Initialize covariance (use process covariance as initial uncertainty)
    m_copy(process_covariance, context->P_km1);

    return context;
}

void kalman_free(kalman_context_t* context) {
    if (!context) return;

    m_free(context->sigma_weights_m);
    m_free(context->sigma_weights_c);

    m_free(context->Q);
    m_free(context->R);

    m_free(context->P_km1);
    m_free(context->sqrt_P_km1);

    m_free(context->P_k_t);
    m_free(context->P_yy);
    m_free(context->P_xy);

    m_free(context->x_hat);
    m_free(context->x_hatm1);
    m_free(context->x_hat_k_t);
    m_free(context->chi_k);
    m_free(context->chi_km1);

    m_free(context->x_b);
    m_free(context->x_a);
    m_free(context->P_a);
    m_free(context->P_b);
    m_free(context->Y_k);
    m_free(context->y_hat_k_t);
    m_free(context->y_a);
    m_free(context->y_b);
    m_free(context->K);
    m_free(context->K_transpose);
    m_free(context->scratch_yx);

    free(context);
}

static error_t _kalman_calc_sigmas(kalman_context_t* context) {
    // First sigma point is the mean
    if (E_OK != m_copy_column(context->x_hat, 0, context->chi_km1, 0)) {
        return E_VAL;
    }

    // Compute Cholesky decomposition of P
    if (E_OK != la_decompositions_cholesky(context->P_km1, context->sqrt_P_km1)) {
        return E_VAL;
    }

    // Generate sigma points
    size_t n_state = context->state_len;
    for (size_t i = 0; i < n_state; i++) {
        // Plus sigma points (indices 1 to n)
        if (E_OK != m_copy_column(context->x_hat, 0, context->chi_km1, i + 1)) {
            return E_VAL;
        }
        if (E_OK != m_add_scaled_column(context->sqrt_P_km1, i, context->eta, context->chi_km1, i + 1)) {
            return E_VAL;
        }

        // Minus sigma points (indices n+1 to 2n)
        if (E_OK != m_copy_column(context->x_hat, 0, context->chi_km1, n_state + i + 1)) {
            return E_VAL;
        }
        if (E_OK != m_add_scaled_column(context->sqrt_P_km1, i, -context->eta, context->chi_km1, n_state + i + 1)) {
            return E_VAL;
        }
    }

    return E_OK;
}

static error_t _kalman_propagate_sigmas(kalman_context_t* context, m_t* input_vector) {
    if (!context || !input_vector) {
        return E_NULLP;
    }

    for (size_t i = 0; i < context->chi_k->cols; i++) {
        // Extract sigma point
        if (E_OK != m_copy_column(context->chi_km1, i, context->x_a, 0)) {
            return E_ERR;
        }

        // Propagate through state function
        if (E_OK != context->state_fn(context->x_a, input_vector, context->x_b)) {
            return E_ERR;
        }

        // Store propagated sigma point
        if (E_OK != m_copy_column(context->x_b, 0, context->chi_k, i)) {
            return E_ERR;
        }
    }

    return E_OK;
}

error_t kalman_step(kalman_context_t* context, m_t* input_vector, m_t* measurement) {
    if (!context || !input_vector || !measurement) {
        return E_NULLP;
    }

    // ============ PREDICTION STEP ============

    // Calculate sigma points from current state estimate
    if (E_OK != _kalman_calc_sigmas(context)) {
        return E_ERR;
    }

    // Propagate sigma points through state transition
    if (E_OK != _kalman_propagate_sigmas(context, input_vector)) {
        return E_ERR;
    }

    // Calculate predicted state mean (x_hat_k_t)
    m_set_all(context->x_hat_k_t, 0);
    for (size_t i = 0; i < context->chi_k->cols; i++) {
        m_data_t w = m_get(context->sigma_weights_m, i, 0);
        m_add_scaled_column(context->chi_k, i, w, context->x_hat_k_t, 0);
    }

    // Calculate predicted state covariance (P_k_t)
    m_set_all(context->P_k_t, 0);
    for (size_t i = 0; i < context->chi_k->cols; i++) {
        // x_a = chi_k[i] - x_hat_k_t
        m_copy_column(context->chi_k, i, context->x_a, 0);
        m_add_scaled_column(context->x_hat_k_t, 0, -1.0, context->x_a, 0);

        // P_a = x_a * x_a^T
        m_outer_product(context->x_a, context->x_a, context->P_a);

        // P_b = weight * P_a
        m_data_t w = m_get(context->sigma_weights_c, i, 0);
        m_scalar_multiply(context->P_a, w, context->P_b);

        // P_k_t += P_b
        m_add(context->P_k_t, context->P_b, context->P_k_t);
    }

    // Add process noise: P_k_t += Q
    m_add(context->P_k_t, context->Q, context->P_k_t);

    // ============ MEASUREMENT UPDATE STEP ============

    // Transform sigma points through measurement function
    for (size_t i = 0; i < context->chi_k->cols; i++) {
        // Extract propagated sigma point
        m_copy_column(context->chi_k, i, context->x_a, 0);

        // Transform through measurement function
        context->measurement_fn(context->x_a, NULL, context->y_a);

        // Store in Y_k
        m_copy_column(context->y_a, 0, context->Y_k, i);
    }

    // Calculate predicted measurement mean (y_hat_k_t)
    m_set_all(context->y_hat_k_t, 0);
    for (size_t i = 0; i < context->Y_k->cols; i++) {
        m_data_t w = m_get(context->sigma_weights_m, i, 0);
        m_add_scaled_column(context->Y_k, i, w, context->y_hat_k_t, 0);
    }

    // Calculate measurement covariance (P_yy)
    m_set_all(context->P_yy, 0);
    for (size_t i = 0; i < context->Y_k->cols; i++) {
        // y_a = Y_k[i] - y_hat_k_t
        m_copy_column(context->Y_k, i, context->y_a, 0);
        m_add_scaled_column(context->y_hat_k_t, 0, -1.0, context->y_a, 0);

        // y_b = y_a * y_a^T (using y_b as scratch for outer product result)
        m_outer_product(context->y_a, context->y_a, context->y_b);

        // Scale by weight
        m_data_t w = m_get(context->sigma_weights_c, i, 0);
        m_scalar_multiply(context->y_b, w, context->y_b);

        // Accumulate
        m_add(context->P_yy, context->y_b, context->P_yy);
    }

    // Add measurement noise: P_yy += R
    m_add(context->P_yy, context->R, context->P_yy);

    // Calculate cross-covariance (P_xy)
    m_set_all(context->P_xy, 0);
    for (size_t i = 0; i < context->chi_k->cols; i++) {
        // x_a = chi_k[i] - x_hat_k_t
        m_copy_column(context->chi_k, i, context->x_a, 0);
        m_add_scaled_column(context->x_hat_k_t, 0, -1.0, context->x_a, 0);

        // y_a = Y_k[i] - y_hat_k_t
        m_copy_column(context->Y_k, i, context->y_a, 0);
        m_add_scaled_column(context->y_hat_k_t, 0, -1.0, context->y_a, 0);

        // Compute weighted outer product and accumulate directly
        for (size_t r = 0; r < context->state_len; r++) {
            for (size_t c = 0; c < context->measurement_len; c++) {
                m_data_t x_val = m_get(context->x_a, r, 0);
                m_data_t y_val = m_get(context->y_a, c, 0);
                m_data_t w = m_get(context->sigma_weights_c, i, 0);
                m_data_t old_val = m_get(context->P_xy, r, c);
                m_set(context->P_xy, r, c, old_val + w * x_val * y_val);
            }
        }
    }

    // Calculate Kalman gain: K = P_xy * P_yy^-1
    // Use the inversion function for positive semi-definite matrices
    // Store P_yy^-1 in y_b (reusing it as measurement x measurement scratch)
    if (E_OK != la_decompositions_invert_positive_semi_definite(context->P_yy, context->y_b)) {
        return E_ERR;
    }

    // K = P_xy * P_yy^-1
    m_mult(context->P_xy, context->y_b, context->K);

    // ============ STATE UPDATE ============

    // innovation = measurement - y_hat_k_t
    // x_hat = x_hat_k_t + K * innovation
    m_copy_column(context->x_hat_k_t, 0, context->x_hat, 0);

    // y_a = measurement - y_hat_k_t
    m_copy_column(measurement, 0, context->y_a, 0);
    m_add_scaled_column(context->y_hat_k_t, 0, -1.0, context->y_a, 0);

    // x_a = K * y_a (innovation)
    m_mult(context->K, context->y_a, context->x_a);

    // x_hat += x_a
    m_add_scaled_column(context->x_a, 0, 1.0, context->x_hat, 0);

    // ============ COVARIANCE UPDATE ============
    // P = P_k_t - K * P_yy * K^T

    m_transpose(context->K, context->K_transpose);

    // scratch_yx = P_yy * K^T
    m_mult(context->P_yy, context->K_transpose, context->scratch_yx);

    // P_a = K * scratch_yx (but P_a is state x state, and K is state x measurement,
    // scratch_yx is measurement x state, so K * scratch_yx is state x state - correct!)
    m_mult(context->K, context->scratch_yx, context->P_a);

    // P_km1 = P_k_t - P_a
    m_scalar_multiply(context->P_a, -1.0, context->P_a);
    m_add(context->P_k_t, context->P_a, context->P_km1);

    return E_OK;
}

error_t kalman_get_state(kalman_context_t* context, m_t* state_out) {
    if (!context || !state_out) {
        return E_NULLP;
    }

    if (state_out->rows != context->state_len || state_out->cols != 1) {
        return E_VAL;
    }

    return m_copy(context->x_hat, state_out);
}

error_t kalman_get_covariance(kalman_context_t* context, m_t* covariance_out) {
    if (!context || !covariance_out) {
        return E_NULLP;
    }

    if (covariance_out->rows != context->state_len || covariance_out->cols != context->state_len) {
        return E_VAL;
    }

    return m_copy(context->P_km1, covariance_out);
}
