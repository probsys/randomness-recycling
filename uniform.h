/*
  Name:     uniform.h
  Purpose:  Generating uniform pseudo-random numbers.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef UNIFORM_H
#define UNIFORM_H

#include "types.h"

struct uniform_preprocessed_s {
    u32 num_outcomes;
    u32 quotient;
    u32 not_remainder;
    u64 inverse;
};

u32 flip(void);
u64 flip_n(u32 n);

extern u64 bits_consumed;
void merge_state(u64 state, u64 bound);
void merge_state_bits(u64 state, u64 n);
u64 uniform_eo(u64 n);
u64 flip_n_from_unif(u32 n);
u32 uniform_u32_from_unif();
bool bernoulli_eo(u32 numer, u32 denom);
struct uniform_preprocessed_s uniform_preprocess(u32 m);
u32 uniform_prediv(struct uniform_preprocessed_s *x);


// Macros for min and max.
#define max(a, b)           \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a > _b ? _a : _b;      \
  })

#define min(a, b)           \
  ({                        \
    __typeof__(a) _a = (a); \
    __typeof__(b) _b = (b); \
    _a < _b ? _a : _b;      \
  })

#endif
