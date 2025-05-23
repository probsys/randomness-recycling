/*
  Name:     aldr.h
  Purpose:  Amplified Loaded Dice Roller.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef ALDR_H
#define ALDR_H

#include "uniform.h"
#include "types.h"

// flattened ALDR tree with entropy-optimal recycling
// (except throwing away accept-reject Bernoulli information)
struct aldr_recycle_s {
    u32 length_breadths;
    u32 length_leaves_flat;
    u32 length_weights;
    u32 reject_weight;
    u32 *breadths;
    u32 *leaves_flat;
    u64 *weights;
};

// FLDR but packing to the left so there is no rejection
struct fldr_eo_s {
  u32 length_breadths;
  u32 length_leaves_flat;
  u32 length_weights;
  struct uniform_preprocessed_s uniform_preprocessed;
  u32 *breadths;
  u32 *leaves_flat;
  u32 *weights;
};

void free_aldr_recycle (struct aldr_recycle_s x);
struct aldr_recycle_s preprocess_aldr_recycle(u32* a, u32 n);
u32 sample_aldr_recycle(struct aldr_recycle_s* f);
u32 bytes_aldr_recycle(struct aldr_recycle_s *x);

void free_fldr_eo(struct fldr_eo_s x);
struct fldr_eo_s preprocess_fldr_eo(u32* a, u32 n);
u32 sample_fldr_eo(struct fldr_eo_s* f);
u32 bytes_fldr_eo(struct fldr_eo_s *x);

#endif
