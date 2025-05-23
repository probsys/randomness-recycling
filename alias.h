/*
  Name:     alias.h
  Purpose:  Exact weighted alias method, translated from Rust to C.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef ALIAS_H
#define ALIAS_H

#include "types.h"

// weighted alias index arrays
struct weighted_alias_s {
    u32 length;
    u32 weight_sum;
    u32 *aliases;
    u32 *no_alias_odds;
};

// weighted alias index arrays with entropy-optimal recycling
struct weighted_alias_eo_s {
    u32 length;
    u32 weight_sum;
    u32 *weights;
    u32 *aliases;
    u32 *no_alias_odds;
    u64 *offsets;
};

void free_weighted_alias(struct weighted_alias_s x);
struct weighted_alias_s preprocess_weighted_alias(int* a, int n);
int bytes_weighted_alias(struct weighted_alias_s *x);

u32 sample_weighted_alias_recycle(struct weighted_alias_s *x);

void free_weighted_alias_eo(struct weighted_alias_eo_s x);
struct weighted_alias_eo_s preprocess_weighted_alias_eo(int* a, int n);
u32 sample_weighted_alias_eo(struct weighted_alias_eo_s *x);
int bytes_weighted_alias_eo(struct weighted_alias_eo_s *x);

#endif
