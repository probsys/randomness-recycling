/*
  Name:     lookup.c
  Purpose:  Table lookup sampling.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>

#include "binarysearch.h"
#include "lookup.h"
#include "types.h"
#include "uniform.h"

struct lookup_eo_s preprocess_lookup_eo(int* a, int n) {
    struct array_s cdf = preprocess_cdf(a, n);
    u32 m = cdf.a[cdf.length - 1];
    struct lookup_eo_s x = {
        .cdf_length = cdf.length,
        .lookup_length = m,
        .cdf = cdf.a,
        .lookup = malloc(m * sizeof(x.lookup[0]))
    };
    for (u32 i = 0; i < n; ++i) {
        for (u32 j = x.cdf[i]; j < x.cdf[i+1]; ++j) {
            x.lookup[j] = i;
        }
    }
    return x;
}

u32 sample_lookup_eo(struct lookup_eo_s *x) {
    u32 uniform_index = uniform_eo(x->lookup_length);
    u32 result = x->lookup[uniform_index];
    merge_state(
        uniform_index - x->cdf[result],
        x->cdf[result + 1] - x->cdf[result]
    );
    return result;
}

void free_lookup_eo(struct lookup_eo_s x) {
    free(x.cdf);
    free(x.lookup);
}

u32 bytes_lookup_eo(struct lookup_eo_s *x) {
    return sizeof(x->cdf_length) +
           sizeof(x->lookup_length) +
           x->cdf_length * sizeof(x->cdf[0]) +
           x->lookup_length * sizeof(x->lookup[0]);
}
