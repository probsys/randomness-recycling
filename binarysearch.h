/*
  Name:     binarysearch.h
  Purpose:  Binary search sampling.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef BINARYSEARCH_H
#define BINARYSEARCH_H

#include "types.h"

struct array_s preprocess_cdf(int* a, int n);
u32 sample_cdf_eo(struct array_s *x);

#endif
