#include <stdint.h>
#include <stdio.h>

/* Includes from the testing source tree */
#include "clar.h"
#include "test.h"

/* Includes from the project source tree */
#include "data_structures/vector.h"


void test_vector__initialize(void) {
    global_test_counter++;
}

void test_vector__cleanup(void) {
}

void test_vector__set(void) {
    v_t *v;
    v = v_new(20);
    v_data_t val = 2.0f;
    cl_assert_(v, "Couldn't create a vector for setting.");

    for (size_t i=0; i < v_len(v); i++)
        cl_assert_equal_i_(v_set(v, i, val), E_OK, "Couldn't set.");

    for (size_t i=0; i < v_len(v); i++)
        cl_assert_(v_get(v,i) == val, "Get returned wrong value after setting.");

    v_del(v);
}
