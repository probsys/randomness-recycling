/*
  Name:     lookup.h
  Purpose:  Table lookup sampling.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#ifndef LOOKUP_H
#define LOOKUP_H

#include "types.h"

struct lookup_eo_s {
    u32 cdf_length;
    u32 lookup_length;
    u32 *cdf;
    u32 *lookup;
};

struct lookup_eo_s preprocess_lookup_eo(int* a, int n);
u32 sample_lookup_eo(struct lookup_eo_s *x);
void free_lookup_eo(struct lookup_eo_s x);
u32 bytes_lookup_eo(struct lookup_eo_s *x);

#endif
