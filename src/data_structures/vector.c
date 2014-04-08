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
