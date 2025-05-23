/*
  Name:     sample.c
  Purpose:  Sampling with randomness recycling.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "uniform.h"
#include "aldr.h"
#include "alias.h"
#include "lookup.h"
#include "binarysearch.h"

#define SAMPLE_PRINT(key, \
        struct_name, \
        func_preprocess, \
        func_sample, \
        func_free) \
    if(strcmp(var_sampler, key) == 0) { \
        struct struct_name s = func_preprocess(a, n); \
        for (u32 i = 0; i < num_samples; ++i) { \
            printf("%d ", func_sample(&s)); \
        } \
        printf("\n"); \
        func_free(s); \
    }

int main(int argc, char **argv) {
    if (argc < 4) {
        printf("usage: %s <sampler (cdf|lookup|alias|fldr|aldr)> <num_samples> <distribution>\n", argv[0]);
        exit(0);
    }
    char *var_sampler = argv[1];
    u32 num_samples = strtoul(argv[2], NULL, 10);

    // Parse the distribution.
    u32 n = argc - 3;
    u32 *a = calloc(n, sizeof(*a));
    for (u32 i = 0; i < n; ++i) {
        a[i] = strtoul(argv[i + 3], NULL, 10);
    }

    // Obtain the samples.
    SAMPLE_PRINT("cdf",
        array_s,
        preprocess_cdf,
        sample_cdf_eo,
        free_array)
    else SAMPLE_PRINT("lookup",
        lookup_eo_s,
        preprocess_lookup_eo,
        sample_lookup_eo,
        free_lookup_eo)
    else SAMPLE_PRINT("alias",
        weighted_alias_eo_s,
        preprocess_weighted_alias_eo,
        sample_weighted_alias_eo,
        free_weighted_alias_eo)
    else SAMPLE_PRINT("fldr",
        fldr_eo_s,
        preprocess_fldr_eo,
        sample_fldr_eo,
        free_fldr_eo)
    else SAMPLE_PRINT("aldr",
        aldr_recycle_s,
        preprocess_aldr_recycle,
        sample_aldr_recycle,
        free_aldr_recycle)

    // Free the heap.
    free(a);

    return 0;
}
