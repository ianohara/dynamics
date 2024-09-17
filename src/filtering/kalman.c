#include "filtering/kalman.h"
#include "linear_algebra/decompositions.h"
#include "data_structures/matrix.h"

kalman_context_t* kalman_new(m_t *initial_state_guess, kalman_state_fn state_fn, m_t* state_covariance, m_t* measurement_covariance) {
    if (!initial_state_guess || !state_fn || !state_covariance || !measurement_covariance) {
        return NULL;
    }

    if (!m_is_vector(initial_state_guess)) {
        return NULL;
    }

    if (!m_is_square(state_covariance) || state_covariance->rows != initial_state_guess->rows) {
        return NULL;
    }

    if (!m_is_square(measurement_covariance)) {
        return NULL;
    }

    size_t state_len = m_max_dim(initial_state_guess);

    kalman_context_t* context = malloc(sizeof(*context));
    // Allocate and zero everything, then calculate our initial cholesky factor
    if (!context) {
        return NULL;
    }

    context->alpha = 1e-4;
    context->beta = 2.0;

    context->lambda = state_len * (context->alpha*context->alpha - 1.0);
    context->eta = sqrtf((state_len + context->lambda));

    context->state_fn = state_fn;
    size_t sigma_point_count = 2*state_len + 1;
    if (!(context->sigma_weights = m_new(sigma_point_count, 1))) {
        kalman_free(context);
        return NULL;
    }

    m_data_t weight_0 = context->lambda / (state_len + context->lambda);
    m_data_t weight_i = 1.0 / (2 * (state_len + context->lambda));
    for (size_t m = 0; m < context->sigma_weights->rows; m++) {
        m_data_t this_weight = m == 0 ? weight_0 : weight_i;
        if (E_OK != m_set(context->sigma_weights, m, 1, this_weight)) {
            kalman_free(context);
            return NULL;
        }
    }

    if (!(context->chi_k = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->chi_km1 = m_new(state_len, sigma_point_count))) {
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

    if (!(context->S_k = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->S_k_t = m_new(state_len, state_len))) {
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

    if (E_OK != la_decompositions_cholesky(initial_covariance, context->S_k)) {
        m_free(initial_covariance);
        kalman_free(context);

        return NULL;
    }


    if (!(context->x_scratch = m_new(state_len, 1))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->weighted_chi_minus_x_k_t = m_new(state_len, sigma_point_count))) {
        kalman_free(context);
        return NULL;
    }

    // Both of the covariance matricies are only used in sqrt form, so just calculate and keep those.
    if (!(context->sqrt_state_covariance = m_new(state_len, state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (E_OK != la_decompositions_cholesky(state_covariance, context->sqrt_state_covariance)) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->sqrt_measurement_covariance = m_new(measurement_covariance->rows, measurement_covariance->cols))) {
        kalman_free(context);
        return NULL;
    }

    if (E_OK != la_decompositions_cholesky(measurement_covariance, context->sqrt_measurement_covariance)) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->state_qr_block = m_new(state_len, 2*state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->state_qr_Q_placeholder = m_new(state_len, 2*state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->measurement_qr_block = m_new(state_len, 2*state_len))) {
        kalman_free(context);
        return NULL;
    }

    if (!(context->measurement_qr_Q_placeholder = m_new(state_len, 2*state_len))) {
        kalman_free(context);
        return NULL;
    }

    return context;
}

void kalman_free(kalman_context_t* context) {
    m_free(context->sigma_weights);
    m_free(context->x_hat);
    m_free(context->x_hatm1);
    m_free(context->x_hat_k_t);
    m_free(context->S_k);
    m_free(context->S_k_t);
    m_free(context->chi_k);
    m_free(context->chi_km1);

    m_free(context->x_scratch);
    m_free(context->weighted_chi_minus_x_k_t);
    m_free(context->sqrt_state_covariance);
    m_free(context->sqrt_measurement_covariance);
    m_free(context->state_qr_block);
    m_free(context->state_qr_Q_placeholder);
    m_free(context->measurement_qr_block);
    m_free(context->measurement_qr_Q_placeholder);
}

static error_t _kalman_calc_sigmas(kalman_context_t* context) {
    if (E_OK != m_copy_column(context->x_hat, 0, context->chi_km1, 0)) {
        return E_VAL;
    }

    for (size_t n = 1; n < context->chi_km1->cols; n++) {
        if (E_OK != m_copy_column(context->x_hat, n, context->chi_km1, n)) {
            return E_VAL;
        }

        if (E_OK != m_add_scaled_column(context->S_k, n, context->eta, context->chi_km1, n)) {
            return E_VAL;
        }

        if (E_OK != m_add_scaled_column(context->S_k, n, -context->eta, context->chi_km1, n)) {
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
        if (E_OK != context->state_fn(context->x_hatm1, input_vector, context->x_scratch)) {
            return E_ERR;
        }

        if (E_OK != m_copy_column(context->x_scratch, 0, context->chi_k, n)) {
            return E_ERR;
        }
    }

    return E_OK;
}

error_t kalman_step(kalman_context_t* context, m_t* input_vector, m_t* measurement) {
    if (!context || !input_vector) {
        return E_NULLP;
    }

    // Equation (17)
    if (E_OK != _kalman_calc_sigmas(context)) {
        return E_ERR;
    }

    // Equation (18)
    if (E_OK != _kalman_propogate_sigmas(context, input_vector)) {
        return E_ERR;
    }

    // Equation (19)
    m_set_all(context->x_hat_k_t, 0);
    for (size_t n = 0; n < context->chi_k->cols; n++) {
        if (E_OK != m_add_scaled_column(context->chi_k, n, m_get(context->sigma_weights, n, 0), context->x_hat_k_t, 0)) {
            return E_ERR;
        }
    }

    // Equation (20)
    // Fill in block matrix for qr decomp
    m_set_all(context->state_qr_block, 0.0);
    m_data_t sqrt_weight_1 = sqrtf(m_get(context->sigma_weights, 1, 0));
    for (size_t n = 0; n < context->chi_k->cols - 1; n++) {
        if (E_OK != m_add_scaled_column(context->chi_k, n+1, sqrt_weight_1, context->state_qr_block, n)) {
            return E_ERR;
        }
        if (E_OK != m_add_scaled_column(context->x_hat_k_t, n+1, -sqrt_weight_1, context->state_qr_block, n)) {
            return E_ERR;
        }
    }

    if (E_OK != m_copy_into(context->sqrt_state_covariance, context->state_qr_block, 0, context->chi_k->cols)) {
        return E_ERR;
    }

    if (E_OK != la_decompositions_qr(context->state_qr_block, context->state_qr_Q_placeholder, context->S_k_t)) {
        return E_ERR;
    }

    // Equation (21)
    m_set_all(context->x_scratch, 0);
    if (E_OK != m_copy_column(context->chi_k, 0, context->x_scratch, 0)) {
        return E_ERR;
    }
    if (E_OK != m_add_scaled_column(context->x_hat_k_t, 0, -1.0, context->x_scratch, 0)) {
        return E_ERR;
    }
    if (E_OK != la_decompositions_cholesky_update(context->S_k_t, context->x_scratch, m_get(context->sigma_weights, 0, 0)) {
        return E_ERR;
    }
}