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

void test_data_structures_vector__cleanup(void)
{
}

void test_data_structures_vector__set(void)
{
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

void test_data_structures_vector__get(void)
{
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

void test_data_structures_vector__dot(void)
{
    size_t len_short = 20, len_long = 1010;
    float flts[len_short], flts2[len_short], flts3[len_short];
    v_t *vfl, *vfl2, *vfl3, *vz, *vo, *vz_long;

    for (size_t i=0; i<array_length(flts); i++) flts[i] = (float)i;
    for (size_t i=0; i<array_length(flts2); i++) flts2[i] = (float)(2*i);
    for (size_t i=0; i<array_length(flts3); i++) flts3[i] = (float)(3*i);

    vfl = v_new_from_floats(flts, array_length(flts));
    vfl2 = v_new_from_floats(flts2, array_length(flts2));
    vfl3 = v_new_from_floats(flts3, array_length(flts3));
    vz = v_new_zeros(len_short);
    vo = v_new_ones(len_short);
    vz_long = v_new_zeros(len_long);

    /* Save vfl3 for later, shut up warnings for now */
    (void)vfl3;

    /*** Input checking ***/
    cl_assert_(v_isnan(v_dot(NULL, NULL)), "v_dot with both vectors NULL should be NAN");
    cl_assert_(v_isnan(v_dot(vz, NULL)), "v_dot with lhs NULL should be NAN");
    cl_assert_(v_isnan(v_dot(NULL, vz)), "v_dot with rhs NULL should be NAN");
    cl_assert_(v_isnan(v_dot(vz_long, vz)), "v_dot with vectors of different length should be NAN");

    /*** Mathematical correctness tests ***/
    cl_assert_(
        v_dot(vfl,vfl2) == v_dot(vfl2,vfl),
        "Commutivity of dot product failed."
        );

    cl_assert_(v_dot(vz, vo) == v_get(vz, 0), "Any dot product with a vector of zeros should be zero!");
    cl_assert_(v_dot(vfl, vz) == v_get(vz, 0), "Any dot product with a vector of zeros should be zero!");
    cl_assert_(v_dot(vo, vo) == (v_data_t)v_len(vo), "Dot product between two vectors of all ones should be length of vector.");

    v_del(vfl);
    v_del(vfl2);
    v_del(vfl3);
    v_del(vz);
    v_del(vo);
    v_del(vz_long);

}
