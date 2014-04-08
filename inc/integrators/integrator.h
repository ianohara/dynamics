#ifndef __INTEGRATOR_H_29948585__
#define __INTEGRATOR_H_29948585__

#include "errors.h"
#include "vector.h"

/* A state_fn is one that takes a state, a control, and calculates the corresponding
   rate of change of state.
*/
typedef error_t (*state_fn)(
                            vector_t *cur_st,
                            vector_t *cur_ctrl,
                            vector_t *cur_st_rate
                            );

/* integrator_fns are those that take a state rate function, a current state, a current control,
   and a next state output pointer and calculate the next state.

   An example integrator_fn would be a 4th order runge-kutta function
*/
typedef error_t (*integrator_fn)(
                                state_fn fn,
                                float    dt,
                                vector_t *cur_st,
                                vector_t *cur_ctrl,
                                vector_t *cur_st_rate
                                );

typedef struct integrator {
    integrator_fn int_fn;
} integrator_t;

integrator_t* integrator_new(integrator_fn int_fn);

#endif /* __INTEGRATOR_H_29948585__ */
