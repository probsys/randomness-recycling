/*
  Name:     types.h
  Purpose:  Shared types.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <stdbool.h>

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

#define u32 uint32_t
#define u64 uint64_t
#define u128 __uint128_t

#define f32 float
#define f64 double

// array
struct array_s
{
    u32 length;
    u32 *a;
};

void free_array(struct array_s x);

u32 bytes_array(struct array_s *x);

#endif