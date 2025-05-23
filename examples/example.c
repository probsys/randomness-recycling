/*
  Name:     example.c
  Purpose:  Example of running ALDR.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>
#include <stdio.h>

#include "aldr.h"
#include "alias.h"
#include "binarysearch.h"
#include "lookup.h"
#include "uniform.h"
#include "types.h"

int main(int argc, char **argv) {
    uint32_t num_samples_each = 18;
    uint32_t num_samplers = 5;
    uint32_t num_samples = num_samples_each * num_samplers;
    uint32_t *samples = calloc(num_samples, sizeof(*samples));

    uint32_t distribution[5] = { 1, 1, 2, 3, 2 };

    struct array_s s_cdf = preprocess_cdf(distribution, 5);
    struct lookup_eo_s s_lookup = preprocess_lookup_eo(distribution, 5);
    struct weighted_alias_eo_s s_alias = preprocess_weighted_alias_eo(distribution, 5);
    struct fldr_eo_s s_fldr = preprocess_fldr_eo(distribution, 5);
    struct aldr_recycle_s s_aldr = preprocess_aldr_recycle(distribution, 5);
    for (uint32_t i = 0; i < num_samples_each; ++i) {
        samples[i] = sample_cdf_eo(&s_cdf);
        samples[i+num_samples_each] = sample_lookup_eo(&s_lookup);
        samples[i+2*num_samples_each] = sample_weighted_alias_eo(&s_alias);
        samples[i+3*num_samples_each] = sample_fldr_eo(&s_fldr);
        samples[i+4*num_samples_each] = sample_aldr_recycle(&s_aldr);
    }
    
    for (u32 i = 0; i < num_samples; ++i) {
        printf("%d ", samples[i]);
    }
    printf("\n");

    free(samples);
    free_array(s_cdf);
    free_lookup_eo(s_lookup);
    free_weighted_alias_eo(s_alias);
    free_fldr_eo(s_fldr);
    free_aldr_recycle(s_aldr);

    return 0;
}
