#include <stdlib.h>

#include "data_structures/vector.h"

v_t* v_new(size_t len) {
    v_t *nv = NULL;
    v_data_t *d = NULL;

    nv = malloc(sizeof *nv);
    if (!nv) goto v_fail;

    d = malloc(len*sizeof *d);
    if (!d) goto v_dat_fail;

    nv->len  = len;
    nv->data = d;
    goto out;

    v_dat_fail:
    free(nv);

    v_fail:
    nv = NULL;

    out:
    return nv;
}

error_t v_del(v_t *v) {
    if (!v) return E_OK;
    if (v->data) free(v->data);

    free(v);
    return E_OK;
}

v_t* v_new_from_floats(float *d, size_t dlen)
{
    v_t *v;
    if (!d || !dlen) return NULL;
    v = v_new(dlen);
    if (!v) return NULL;

    for (size_t i=0; i < dlen; i++)
        v_set(v, i, (v_data_t)d[i]);

    return v;
}

v_t* v_new_from_value(v_data_t val, size_t len)
{
    v_t *v = NULL;
    if (!len) return NULL;

    v = v_new(len);
    if (!v) return NULL;

    for (size_t i=0; i < v_len(v); i++)
        v_set(v,i,val);

    return v;
}

v_t* v_new_zeros(size_t len) {
    return v_new_from_value((v_data_t)0.0, len);
}

v_t* v_new_ones(size_t len) {
    return v_new_from_value((v_data_t)1.0, len);
}

error_t v_set(const v_t *v, size_t ind, v_data_t val) {
    if (!v) return E_NULLP;
    if (ind >= v->len) return E_VAL;
    v->data[ind] = val;
    return E_OK;
}

v_data_t v_get(const v_t * const v, size_t ind) {
    if (!v) return V_NAN;
    if (ind >= v->len) return V_NAN;
    return v->data[ind];
}

size_t v_len(const v_t * const v) {
    if (!v) return 0;
    return v->len;
}

error_t v_sp(v_data_t s, v_t *v, v_t *res)
{
    if (!v || !res) return E_NULLP;
    if (v_len(v) != v_len(res)) return E_VAL;

    for (size_t i=0; i < v_len(res); i++)
        v_set(res,i,s*v_get(v,i));

    return E_OK;
}

error_t v_sum(v_t *lhs, v_t *rhs, v_t *res) {
    if (!lhs || !rhs || !res) return E_NULLP;
    if (lhs->len != rhs->len) return E_VAL;
    if (res->len != lhs->len) return E_VAL;

    for (size_t i=0; i < res->len; i++)
        v_set(res, i, v_get(lhs,i)+v_get(rhs,i));

    return E_OK;
}

error_t v_negate(v_t *v, v_t *res) {
    if (!v || !res) return E_NULLP;
    if (v->len != res->len) return E_VAL;

    for (size_t i=0; i < res->len; i++)
        v_set(res, i, -v_get(v,i));

    return E_OK;
}

v_data_t v_dot(v_t *lhs, v_t *rhs) {
    v_data_t res = 0.0;
    if (!lhs || !rhs) return V_NAN;
    if (v_len(lhs) != v_len(rhs)) return V_NAN;

    for (size_t i=0; i < v_len(lhs); i++)
        res += v_get(lhs,i)*v_get(rhs,i);

    return res;
}
