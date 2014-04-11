#include <stdint.h>
#include <stdio.h>
#include <math.h>

/* Includes from the testing source tree */
#include "clar.h"
#include "test.h"

/* Includes from the project source tree */
#include "data_structures/vector.h"

void test_data_structures_vector__initialize(void) {
    global_test_counter++;
}

void test_data_structures_vector__cleanup(void) {
}

void test_data_structures_vector__set(void) {
    v_t *v;
    v = v_new(20);
    v_data_t val = 2.0f;
    cl_assert_(v, "Couldn't create a vector for setting.");

    cl_assert_(v_set(NULL,1,0.0f) == E_NULLP, "v_set on null pointer should fail.");
    cl_assert_(v_set(v, v_len(v)+1,0.0f) == E_VAL, "v_set on index > length shoudl fail.");

    for (size_t i=0; i < v_len(v); i++)
        cl_assert_equal_i_(v_set(v, i, val), E_OK, "Couldn't set.");

    for (size_t i=0; i < v_len(v); i++)
        cl_assert_(v_get(v,i) == val, "Get returned wrong value after setting.");

    v_del(v);
}

void test_data_structures_vector__get(void) {
    v_t *v;
    v = v_new(20);
    cl_assert(v);

    cl_assert_(isnan(v_get(NULL,0)), "v_get on null pointer should fail with NAN");
    cl_assert_(isnan(v_get(v, v_len(v)+1)), "v_get with index that is too high should fail.");

    for (size_t i=0; i < v_len(v); i++)
        cl_assert_equal_i_(v_set(v,i,(v_data_t)i), E_OK, "Couldn't set.");

    for (size_t i=0; i<v_len(v); i++)
        cl_assert_(v_get(v,i) == (v_data_t)i, "Getting failed with wrong value!");

}
