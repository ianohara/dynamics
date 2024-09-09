#ifndef __FILTERING_KALMAN_H__
#define __FILTERING_KALMAN_H__

#include "errors.h"
#include "data_structures/matrix.h"

typedef struct kalman_context {
    /* Alpha must hold to: 1e-4 <= alpha <= 1*/
    double alpha = 1e-4;

    /* Beta of 2 is best for gaussian distributed state vectors */
    double beta = 2;

    // Any of the fields with leading underscores are internal scratch pad values that you
    // should not touch.
    double _lamda;
    double _eta;
    m_t* _sigma_weights;

} kalman_context_t;

#endif