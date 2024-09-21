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

    size_t state_len = m_max_dim(initial_state_guess) + m_max_dim(process_covariance) + m_max_dim(measurement_covariance);
    size_t measurement_len = m_max_dim(measurement_covariance);

    kalman_context_t* context = malloc(sizeof(*context));
    // Allocate and zero everything, then calculate our initial cholesky factor
    if (!context) {
        return NULL;
    }

    context->alpha = 1e-4;
    context->beta = 2.0;

    context->lambda = context->alpha*context->alpha * (state_len + context->kappa) - state_len;
    context->eta = sqrtf(state_len + context->lambda);

    context->state_fn = state_fn;
    context->measurement_fn = measurement_fn;
    size_t sigma_point_count = 2*state_len + 1;
    if (!(context->sigma_weights_m = m_new(sigma_point_count, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->sigma_weights_c = m_new(sigma_point_count, 1))) {
        kalman_free(context);
        return NULL;
    }

    m_data_t weight_0_m = context->lambda / (state_len + context->lambda);
    m_data_t weight_0_c = context->lambda / (state_len + context->lambda) + (1 - context->alpha*context->alpha + context->beta);
    m_data_t weight_i = 1.0 / (2.0 * (state_len + context->lambda));
    for (size_t m = 0; m < context->sigma_weights_m->rows; m++) {
        m_data_t this_weight_m = m == 0 ? weight_0_m : weight_i;
        m_data_t this_weight_c = m == 0 ? weight_0_c : weight_i;
        if (E_OK != m_set(context->sigma_weights_m, m, 1, this_weight_m)) {
            kalman_free(context);
            return NULL;
        }

        if (E_OK != m_set(context->sigma_weights_c, m, 1, this_weight_c)) {
            kalman_free(context);
            return NULL;
        }
    }


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

    if (!(context->chi_k = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->chi_km1 = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

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

    m_t* guess_negated = m_new(state_len, 1);
    if (E_OK != m_negate(initial_state_guess, guess_negated)) {
        m_free(guess_negated);
        kalman_free(context);
        return NULL;
    }

    m_t* initial_covariance = m_new(state_len, state_len);
    if (E_OK != m_outer_product(guess_negated, guess_negated, initial_covariance)) {
        m_free(guess_negated);
        m_free(initial_covariance);
        kalman_free(context);
        return NULL;
    }
    m_free(guess_negated);

    return context;
}

void kalman_free(kalman_context_t* context) {
    m_free(context->sigma_weights_m);
    m_free(context->sigma_weights_c);

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
    m_free(context->Y_k);
    m_free(context->y_hat_k_t);
    m_free(context->y_a);
    m_free(context->K);
    m_free(context->K_transpose);
    m_free(context->scratch_yx);

    free(context);
}

static error_t _kalman_calc_sigmas(kalman_context_t* context) {
    if (E_OK != m_copy_column(context->x_hat, 0, context->chi_km1, 0)) {
        return E_VAL;
    }

    if (E_OK != la_decompositions_cholesky(context->P_km1, context->sqrt_P_km1)) {
        return E_VAL;
    }

    for (size_t n = 1; n < (context->chi_km1->cols - 1) / 2; n++) {
        if (E_OK != m_copy_column(context->x_hat, n, context->chi_km1, n)) {
            return E_VAL;
        }

        if (E_OK != m_add_scaled_column(context->sqrt_P_km1, n, context->eta, context->chi_km1, n)) {
            return E_VAL;
        }
    }

    for (size_t n = (context->chi_km1->cols - 1); n < context->chi_km1->cols; n++) {
        if (E_OK != m_copy_column(context->x_hat, n, context->chi_km1, n)) {
            return E_VAL;
        }

        if (E_OK != m_add_scaled_column(context->sqrt_P_km1, n, -context->eta, context->chi_km1, n)) {
            return E_VAL;
        }
    }

    return E_OK;
}

static error_t _kalman_propogate_sigmas(kalman_context_t* context, m_t* input_vector) {
    if (!context || !input_vector) {
        return E_NULLP;
    }

    for (size_t n = 0; n < context->chi_k->cols; n++) {
        if (E_OK != m_copy_column(context->chi_km1, n, context->x_a, 0)) {
            return E_ERR;
        }

        // This is wrong --> need to unpack x_a into the x and v parts, and figure out how
        // input vector interplays with v part (they are the same?)
        if (E_OK != context->state_fn(context->x_a, input_vector, context->x_b)) {
            return E_ERR;
        }

        if (E_OK != m_copy_column(context->x_b, 0, context->chi_k, n)) {
            return E_ERR;
        }
    }

    return E_OK;
}

error_t kalman_step(kalman_context_t* context, m_t* input_vector, m_t* measurement) {
    if (!context || !input_vector) {
        return E_NULLP;
    }

    // Calculate chi_k minus 1
    if (E_OK != _kalman_calc_sigmas(context)) {
        return E_ERR;
    }

    //// Time update steps
    // Calculate chi_k
    if (E_OK != _kalman_propogate_sigmas(context, input_vector)) {
        return E_ERR;
    }

    // Calculate x_hat_k_bar (aka x_hat_k_t here)
    m_set_all(context->x_hat_k_t, 0);
    for (size_t n = 0; n < context->chi_k->cols; n++) {
        if (E_OK != m_add_scaled_column(context->chi_k, n, m_get(context->sigma_weights_m, n, 0), context->x_hat_k_t, 0)) {
            return E_ERR;
        }
    }

    // Calculate the time covariance
    m_set_all(context->P_k_t, 0);
    for (size_t n = 0; n < context->chi_k->cols; n++) {
        m_copy_column(context->chi_k, n, context->x_a, 0);
        m_add_scaled_column(context->x_hat_k_t, 0, -1.0, context->x_a, 0);

        m_outer_product(context->x_a, context->x_a, context->P_a);
        m_scalar_multiply(context->P_a, m_get(context->sigma_weights_c, n, 0), context->P_a);

        m_add(context->P_k_t, context->P_a, context->P_k_t);
    }

    // Calculate the measurement estimate given these sigmas
    for (size_t n = 0; n < context->Y_k->cols; n++) {

        context->measurement_fn(context->x_hat_k_t, NULL, context->y_a);
        m_copy_column(context->y_a, 0, context->Y_k, n);
    }


    m_set_all(context->y_hat_k_t, 0);
    for (size_t n = 0; n < context->Y_k->cols; n++) {
        m_add_scaled_column(context->Y_k, n, m_get(context->sigma_weights_m, n, 0), context->y_hat_k_t, 0);
    }

    //// Measurement update calcs
    // covariance of measurements from sigma points
    m_set_all(context->P_yy, 0);
    for (size_t n = 0; n < context->Y_k->cols; n++) {
        m_copy_column(context->Y_k, n, context->y_a, 0);
        m_add_scaled_column(context->y_hat_k_t, 0, -1.0, context->y_a, 0);

        m_outer_product(context->y_a, context->y_a, context->P_yy);
        m_scalar_multiply(context->P_yy, m_get(context->sigma_weights_c, n, 0), context->P_yy);
    }

    // Covariance of the state and sigma estimated measurements
    // NOTE(imo): Could re-use x_a and y_a from above, but this is already confusing so just redo to stay 1-1 w the paper
    m_set_all(context->P_xy, 0);
    for (size_t n = 0; n < context->chi_k->cols; n++) {
        m_copy_column(context->chi_k, n, context->x_a, 0);
        m_add_scaled_column(context->x_hat_k_t, 0, -1.0, context->x_a, 0);

        m_copy_column(context->Y_k, n, context->y_a, 0);
        m_add_scaled_column(context->y_hat_k_t, 0, -1.0, context->y_a, 0);

        m_outer_product(context->x_a, context->y_a, context->P_xy);
        m_scalar_multiply(context->P_xy, m_get(context->sigma_weights_c, n, 0), context->P_xy);
    }

    m_invert(context->P_yy, context->P_a);

    m_mult(context->P_xy, context->P_a, context->K);

    m_copy_column(context->x_hat_k_t, 0, context->x_hat, 0);

    m_mult(context->K, measurement, context->x_a);
    m_mult(context->K, context->y_hat_k_t, context->x_b);
    m_add_scaled_column(context->x_a, 0, 1.0, context->x_hat, 0);
    m_add_scaled_column(context->x_b, 0, -1.0, context->x_hat, 0);

    m_transpose(context->K, context->K_transpose);
    m_mult(context->P_yy, context->K_transpose, context->scratch_yx);

    m_mult(context->K, context->scratch_yx, context->P_a);
    m_scalar_multiply(context->P_a, -1.0, context->P_a);

    // Really P_k in the paper, but just store it immediately in P_km1 since that's how it'll
    // be used next time around.
    m_add(context->P_k_t, context->P_a, context->P_km1);

    return E_OK;
}