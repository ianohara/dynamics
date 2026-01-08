#include <string.h>

#include "filtering/kalman.h"
#include "linear_algebra/decompositions.h"
#include "data_structures/matrix.h"

kalman_context_t* kalman_new(
    m_t *initial_state_guess,
    m_t *initial_covariance,
    kalman_state_fn state_fn,
    kalman_state_to_measurement_fn measurement_fn,
    m_t* process_covariance,
    m_t* measurement_covariance) {
    if (!initial_state_guess || !state_fn || !process_covariance || !measurement_covariance) {
        return NULL;
    }

    // If initial_covariance provided, validate it
    if (initial_covariance) {
        if (!m_is_square(initial_covariance) || initial_covariance->rows != initial_state_guess->rows) {
            return NULL;
        }
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

    // Initialize covariance
    if (initial_covariance) {
        m_copy(initial_covariance, context->P_km1);
    } else {
        // Default to process covariance if no initial covariance given
        m_copy(process_covariance, context->P_km1);
    }

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

// ============================================================================
// Square-Root Unscented Kalman Filter Implementation
// ============================================================================

kalman_sqrt_context_t* kalman_sqrt_new(
    m_t *initial_state_guess,
    m_t *initial_covariance,
    kalman_state_fn state_fn,
    kalman_state_to_measurement_fn measurement_fn,
    m_t *process_covariance,
    m_t *measurement_covariance) {

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

    if (initial_covariance) {
        if (!m_is_square(initial_covariance) || initial_covariance->rows != initial_state_guess->rows) {
            return NULL;
        }
    }

    size_t n = initial_state_guess->rows;
    size_t m = measurement_covariance->rows;
    size_t num_sigma = 2 * n + 1;

    kalman_sqrt_context_t* ctx = malloc(sizeof(*ctx));
    if (!ctx) {
        return NULL;
    }
    memset(ctx, 0, sizeof(*ctx));

    ctx->state_len = n;
    ctx->measurement_len = m;

    // UKF parameters
    ctx->alpha = 1e-3;
    ctx->beta = 2.0;
    ctx->kappa = 0.0;

    ctx->lambda = ctx->alpha * ctx->alpha * (n + ctx->kappa) - n;
    ctx->gamma = sqrt(n + ctx->lambda);

    ctx->state_fn = state_fn;
    ctx->measurement_fn = measurement_fn;

    // Allocate weights
    ctx->Wm = m_new(num_sigma, 1);
    ctx->Wc = m_new(num_sigma, 1);
    if (!ctx->Wm || !ctx->Wc) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Calculate weights
    m_data_t Wm_0 = ctx->lambda / (n + ctx->lambda);
    ctx->Wc_0 = Wm_0 + (1 - ctx->alpha * ctx->alpha + ctx->beta);
    m_data_t Wi = 1.0 / (2.0 * (n + ctx->lambda));

    m_set(ctx->Wm, 0, 0, Wm_0);
    m_set(ctx->Wc, 0, 0, ctx->Wc_0);
    for (size_t i = 1; i < num_sigma; i++) {
        m_set(ctx->Wm, i, 0, Wi);
        m_set(ctx->Wc, i, 0, Wi);
    }

    // Allocate sqrt of noise covariances
    ctx->sqrt_Q = m_new(n, n);
    ctx->sqrt_R = m_new(m, m);
    if (!ctx->sqrt_Q || !ctx->sqrt_R) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Compute sqrt(Q) and sqrt(R) via Cholesky
    if (E_OK != la_decompositions_cholesky(process_covariance, ctx->sqrt_Q)) {
        kalman_sqrt_free(ctx);
        return NULL;
    }
    if (E_OK != la_decompositions_cholesky(measurement_covariance, ctx->sqrt_R)) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Allocate state estimate and sqrt covariance
    ctx->x_hat = m_new(n, 1);
    ctx->S = m_new(n, n);
    if (!ctx->x_hat || !ctx->S) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Initialize state
    m_copy(initial_state_guess, ctx->x_hat);

    // Initialize S (sqrt of covariance)
    m_t* init_cov = initial_covariance ? initial_covariance : process_covariance;
    if (E_OK != la_decompositions_cholesky(init_cov, ctx->S)) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Allocate sigma points
    ctx->chi = m_new(n, num_sigma);
    ctx->chi_prop = m_new(n, num_sigma);
    if (!ctx->chi || !ctx->chi_prop) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Allocate measurement matrices
    ctx->Y = m_new(m, num_sigma);
    ctx->y_hat = m_new(m, 1);
    if (!ctx->Y || !ctx->y_hat) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Allocate cross-covariance and gain
    ctx->Pxy = m_new(n, m);
    ctx->Sy = m_new(m, m);
    ctx->K = m_new(n, m);
    if (!ctx->Pxy || !ctx->Sy || !ctx->K) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // Allocate scratch matrices
    ctx->x_scratch = m_new(n, 1);
    ctx->y_scratch = m_new(m, 1);
    if (!ctx->x_scratch || !ctx->y_scratch) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    // QR matrices - need enough rows for compound matrix
    // For time update: (2n) x n matrix
    // For measurement update: (2n + m) x m matrix (larger)
    size_t qr_rows = 2 * n + m;
    size_t qr_cols = (n > m) ? n : m;
    ctx->qr_input = m_new(qr_rows, qr_cols);
    ctx->qr_Q = m_new(qr_rows, qr_rows);
    ctx->qr_R = m_new(qr_rows, qr_cols);
    if (!ctx->qr_input || !ctx->qr_Q || !ctx->qr_R) {
        kalman_sqrt_free(ctx);
        return NULL;
    }

    return ctx;
}

void kalman_sqrt_free(kalman_sqrt_context_t *ctx) {
    if (!ctx) return;

    m_free(ctx->Wm);
    m_free(ctx->Wc);
    m_free(ctx->sqrt_Q);
    m_free(ctx->sqrt_R);
    m_free(ctx->x_hat);
    m_free(ctx->S);
    m_free(ctx->chi);
    m_free(ctx->chi_prop);
    m_free(ctx->Y);
    m_free(ctx->y_hat);
    m_free(ctx->Pxy);
    m_free(ctx->Sy);
    m_free(ctx->K);
    m_free(ctx->x_scratch);
    m_free(ctx->y_scratch);
    m_free(ctx->qr_input);
    m_free(ctx->qr_Q);
    m_free(ctx->qr_R);

    free(ctx);
}

// Helper: Generate sigma points from x and S
static error_t _sqrt_generate_sigma_points(kalman_sqrt_context_t* ctx) {
    size_t n = ctx->state_len;

    // First sigma point is the mean
    m_copy_column(ctx->x_hat, 0, ctx->chi, 0);

    // Sigma points 1 to n: x + gamma * S[:,i]
    for (size_t i = 0; i < n; i++) {
        m_copy_column(ctx->x_hat, 0, ctx->chi, i + 1);
        m_add_scaled_column(ctx->S, i, ctx->gamma, ctx->chi, i + 1);
    }

    // Sigma points n+1 to 2n: x - gamma * S[:,i]
    for (size_t i = 0; i < n; i++) {
        m_copy_column(ctx->x_hat, 0, ctx->chi, n + i + 1);
        m_add_scaled_column(ctx->S, i, -ctx->gamma, ctx->chi, n + i + 1);
    }

    return E_OK;
}

// Helper: Compute weighted mean of columns
static error_t _sqrt_weighted_mean(m_t* points, m_t* weights, m_t* mean) {
    m_set_all(mean, 0);
    for (size_t i = 0; i < points->cols; i++) {
        m_data_t w = m_get(weights, i, 0);
        m_add_scaled_column(points, i, w, mean, 0);
    }
    return E_OK;
}

error_t kalman_sqrt_step(kalman_sqrt_context_t *ctx, m_t *input_vector, m_t *measurement) {
    if (!ctx || !input_vector || !measurement) {
        return E_NULLP;
    }

    size_t n = ctx->state_len;
    size_t m = ctx->measurement_len;
    size_t num_sigma = 2 * n + 1;

    // ========== PREDICTION STEP ==========

    // Generate sigma points
    _sqrt_generate_sigma_points(ctx);

    // Propagate sigma points through state function
    m_t* prop_out = m_new(n, 1);
    if (!prop_out) return E_ERR;

    for (size_t i = 0; i < num_sigma; i++) {
        m_copy_column(ctx->chi, i, ctx->x_scratch, 0);
        ctx->state_fn(ctx->x_scratch, input_vector, prop_out);
        m_copy_column(prop_out, 0, ctx->chi_prop, i);
    }
    m_free(prop_out);

    // Compute predicted mean
    _sqrt_weighted_mean(ctx->chi_prop, ctx->Wm, ctx->x_hat);

    // Compute predicted sqrt covariance using QR decomposition
    // Form compound matrix: [ sqrt(Wi) * (chi_prop[:,i] - x_hat)^T ]
    //                       [ sqrt_Q^T                             ]
    // Then do QR, and S_pred = R^T (the R from QR of this compound matrix transposed)

    // For simplicity, use a different approach: direct Cholesky updates
    // Start with sqrt_Q, then do rank-1 updates for each sigma point deviation

    // Actually, let's use a simpler compound matrix approach
    // Build matrix where columns are sqrt(Wi)*(chi_i - x_hat) for i=1..2n, plus sqrt_Q columns
    // Then take QR and extract upper triangular R

    // Simpler approach: compute predicted P, then Cholesky
    // P_pred = sum Wi * (chi_i - x_hat)(chi_i - x_hat)^T + Q

    m_t* P_pred = m_new(n, n);
    m_t* outer = m_new(n, n);
    if (!P_pred || !outer) {
        m_free(P_pred);
        m_free(outer);
        return E_ERR;
    }

    m_set_all(P_pred, 0);
    for (size_t i = 0; i < num_sigma; i++) {
        m_copy_column(ctx->chi_prop, i, ctx->x_scratch, 0);
        m_add_scaled_column(ctx->x_hat, 0, -1.0, ctx->x_scratch, 0);

        m_outer_product(ctx->x_scratch, ctx->x_scratch, outer);
        m_data_t w = m_get(ctx->Wc, i, 0);
        m_scalar_multiply(outer, w, outer);
        m_add(P_pred, outer, P_pred);
    }

    // Add Q
    m_t* Q_full = m_new(n, n);
    if (!Q_full) {
        m_free(P_pred);
        m_free(outer);
        return E_ERR;
    }
    // Q = sqrt_Q * sqrt_Q^T
    m_t* sqrt_Q_T = m_new(n, n);
    if (!sqrt_Q_T) {
        m_free(P_pred);
        m_free(outer);
        m_free(Q_full);
        return E_ERR;
    }
    m_transpose(ctx->sqrt_Q, sqrt_Q_T);
    m_mult(ctx->sqrt_Q, sqrt_Q_T, Q_full);
    m_add(P_pred, Q_full, P_pred);
    m_free(sqrt_Q_T);
    m_free(Q_full);

    // Compute sqrt of P_pred
    if (E_OK != la_decompositions_cholesky(P_pred, ctx->S)) {
        m_free(P_pred);
        m_free(outer);
        return E_ERR;
    }
    m_free(P_pred);

    // ========== MEASUREMENT UPDATE STEP ==========

    // Transform sigma points through measurement function
    for (size_t i = 0; i < num_sigma; i++) {
        m_copy_column(ctx->chi_prop, i, ctx->x_scratch, 0);
        ctx->measurement_fn(ctx->x_scratch, NULL, ctx->y_scratch);
        m_copy_column(ctx->y_scratch, 0, ctx->Y, i);
    }

    // Compute predicted measurement mean
    _sqrt_weighted_mean(ctx->Y, ctx->Wm, ctx->y_hat);

    // Compute Pyy (measurement covariance) and take sqrt
    m_t* Pyy = m_new(m, m);
    m_t* outer_y = m_new(m, m);
    if (!Pyy || !outer_y) {
        m_free(Pyy);
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }

    m_set_all(Pyy, 0);
    for (size_t i = 0; i < num_sigma; i++) {
        m_copy_column(ctx->Y, i, ctx->y_scratch, 0);
        m_add_scaled_column(ctx->y_hat, 0, -1.0, ctx->y_scratch, 0);

        m_outer_product(ctx->y_scratch, ctx->y_scratch, outer_y);
        m_data_t w = m_get(ctx->Wc, i, 0);
        m_scalar_multiply(outer_y, w, outer_y);
        m_add(Pyy, outer_y, Pyy);
    }

    // Add R
    m_t* R_full = m_new(m, m);
    m_t* sqrt_R_T = m_new(m, m);
    if (!R_full || !sqrt_R_T) {
        m_free(Pyy);
        m_free(outer_y);
        m_free(outer);
        m_free(R_full);
        m_free(sqrt_R_T);
        return E_ERR;
    }
    m_transpose(ctx->sqrt_R, sqrt_R_T);
    m_mult(ctx->sqrt_R, sqrt_R_T, R_full);
    m_add(Pyy, R_full, Pyy);
    m_free(sqrt_R_T);
    m_free(R_full);

    // Sy = chol(Pyy)
    if (E_OK != la_decompositions_cholesky(Pyy, ctx->Sy)) {
        m_free(Pyy);
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }
    m_free(Pyy);

    // Compute cross-covariance Pxy
    m_set_all(ctx->Pxy, 0);
    for (size_t i = 0; i < num_sigma; i++) {
        m_copy_column(ctx->chi_prop, i, ctx->x_scratch, 0);
        m_add_scaled_column(ctx->x_hat, 0, -1.0, ctx->x_scratch, 0);

        m_copy_column(ctx->Y, i, ctx->y_scratch, 0);
        m_add_scaled_column(ctx->y_hat, 0, -1.0, ctx->y_scratch, 0);

        // Add w * x_scratch * y_scratch^T to Pxy
        m_data_t w = m_get(ctx->Wc, i, 0);
        for (size_t r = 0; r < n; r++) {
            for (size_t c = 0; c < m; c++) {
                m_data_t val = m_get(ctx->Pxy, r, c);
                val += w * m_get(ctx->x_scratch, r, 0) * m_get(ctx->y_scratch, c, 0);
                m_set(ctx->Pxy, r, c, val);
            }
        }
    }

    // Compute Kalman gain: K = Pxy * Sy^-T * Sy^-1 = Pxy * (Sy * Sy^T)^-1
    // Or equivalently, solve Sy * Sy^T * K^T = Pxy^T for K

    // For simplicity, compute Pyy^-1 and multiply
    m_t* Sy_inv = m_new(m, m);
    if (!Sy_inv) {
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }

    // Invert Sy (lower triangular)
    m_set_all(Sy_inv, 0);
    for (size_t j = 0; j < m; j++) {
        for (size_t i = 0; i < m; i++) {
            if (i < j) {
                m_set(Sy_inv, i, j, 0);
            } else if (i == j) {
                m_set(Sy_inv, i, j, 1.0 / m_get(ctx->Sy, i, i));
            } else {
                m_data_t sum = 0;
                for (size_t k = j; k < i; k++) {
                    sum += m_get(ctx->Sy, i, k) * m_get(Sy_inv, k, j);
                }
                m_set(Sy_inv, i, j, -sum / m_get(ctx->Sy, i, i));
            }
        }
    }

    // Pyy_inv = Sy_inv^T * Sy_inv
    m_t* Sy_inv_T = m_new(m, m);
    m_t* Pyy_inv = m_new(m, m);
    if (!Sy_inv_T || !Pyy_inv) {
        m_free(Sy_inv);
        m_free(Sy_inv_T);
        m_free(Pyy_inv);
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }
    m_transpose(Sy_inv, Sy_inv_T);
    m_mult(Sy_inv_T, Sy_inv, Pyy_inv);

    // K = Pxy * Pyy_inv
    m_mult(ctx->Pxy, Pyy_inv, ctx->K);

    m_free(Sy_inv);
    m_free(Sy_inv_T);
    m_free(Pyy_inv);

    // ========== STATE UPDATE ==========

    // innovation = measurement - y_hat
    m_copy(measurement, ctx->y_scratch);
    m_add_scaled_column(ctx->y_hat, 0, -1.0, ctx->y_scratch, 0);

    // x_hat = x_hat + K * innovation
    m_mult(ctx->K, ctx->y_scratch, ctx->x_scratch);
    m_add_scaled_column(ctx->x_scratch, 0, 1.0, ctx->x_hat, 0);

    // ========== COVARIANCE UPDATE ==========
    // S_new = cholupdate(S, K*Sy, '-') for each column of K*Sy

    m_t* U = m_new(n, m);
    if (!U) {
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }
    m_mult(ctx->K, ctx->Sy, U);

    // Downdate S with each column of U
    m_t* u_col = m_new(n, 1);
    if (!u_col) {
        m_free(U);
        m_free(outer_y);
        m_free(outer);
        return E_ERR;
    }

    for (size_t j = 0; j < m; j++) {
        m_copy_column(U, j, u_col, 0);
        if (E_OK != la_decompositions_cholesky_downdate(ctx->S, u_col)) {
            // Downdate failed - fall back to recomputing S from scratch
            // Compute P = S*S^T - U*U^T, then S = chol(P)
            m_t* S_T = m_new(n, n);
            m_t* P_new = m_new(n, n);
            m_t* U_T = m_new(m, n);
            m_t* UUT = m_new(n, n);
            if (S_T && P_new && U_T && UUT) {
                m_transpose(ctx->S, S_T);
                m_mult(ctx->S, S_T, P_new);
                m_transpose(U, U_T);
                m_mult(U, U_T, UUT);
                m_scalar_multiply(UUT, -1.0, UUT);
                m_add(P_new, UUT, P_new);
                la_decompositions_cholesky(P_new, ctx->S);
            }
            m_free(S_T);
            m_free(P_new);
            m_free(U_T);
            m_free(UUT);
            break;
        }
    }

    m_free(u_col);
    m_free(U);
    m_free(outer_y);
    m_free(outer);

    return E_OK;
}

error_t kalman_sqrt_get_state(kalman_sqrt_context_t *ctx, m_t *state_out) {
    if (!ctx || !state_out) {
        return E_NULLP;
    }
    if (state_out->rows != ctx->state_len || state_out->cols != 1) {
        return E_VAL;
    }
    return m_copy(ctx->x_hat, state_out);
}

error_t kalman_sqrt_get_covariance(kalman_sqrt_context_t *ctx, m_t *covariance_out) {
    if (!ctx || !covariance_out) {
        return E_NULLP;
    }
    if (covariance_out->rows != ctx->state_len || covariance_out->cols != ctx->state_len) {
        return E_VAL;
    }

    // P = S * S^T
    m_t* S_T = m_new(ctx->state_len, ctx->state_len);
    if (!S_T) {
        return E_ERR;
    }
    m_transpose(ctx->S, S_T);
    m_mult(ctx->S, S_T, covariance_out);
    m_free(S_T);
    return E_OK;
}

error_t kalman_sqrt_get_sqrt_covariance(kalman_sqrt_context_t *ctx, m_t *sqrt_cov_out) {
    if (!ctx || !sqrt_cov_out) {
        return E_NULLP;
    }
    if (sqrt_cov_out->rows != ctx->state_len || sqrt_cov_out->cols != ctx->state_len) {
        return E_VAL;
    }
    return m_copy(ctx->S, sqrt_cov_out);
}
