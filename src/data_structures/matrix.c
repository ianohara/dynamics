#include <stdlib.h>

#include "data_structures/matrix.h"

m_t* m_new(size_t rows, size_t cols) {
    m_t *nm = NULL;
    m_data_t *d = NULL;
    if (rows <= 0) goto fail;
    if (cols <= 0) goto fail;

    nm = malloc(sizeof *nm);
    if (!nm) goto fail_nm;

    d = malloc(rows*cols*sizeof *d);
    if (!d) goto fail_d;

    nm->rows = rows;
    nm->cols = cols;
    nm->data = d;
    goto out;

    fail_d:
    free(nm);

    fail_nm:
    fail:
    nm = NULL;

    out:
    return nm;
}

error_t m_del(m_t *m) {
    if (!m) return E_OK;
    if (m->data) free(m->data);
    free(m);
    return E_OK;
}
