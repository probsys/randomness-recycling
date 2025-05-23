/*
  Name:     binarysearch.c
  Purpose:  Binary search sampling.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>

#include "binarysearch.h"
#include "types.h"
#include "uniform.h"

struct array_s preprocess_cdf(int* a, int n) {
    struct array_s x = { .length = n+1, .a = malloc((n + 1) * sizeof(u32)) };
    x.a[0] = 0;
    for (u32 i = 0; i < n; ++i) {
        x.a[i + 1] = x.a[i] + a[i];
    }
    return x;
}

u32 sample_cdf_eo(struct array_s *x) {
    u32 uniform_index = uniform_eo(x->a[x->length - 1]);
    u32 low = 1;
    u32 high = x->length - 1;
    while (low < high) {
        u32 mid = (low + high) / 2;
        if (x->a[mid] <= uniform_index) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    merge_state(uniform_index - x->a[low-1], x->a[low] - x->a[low-1]);
    return low - 1;
}
