/*
  Name:     types.c
  Purpose:  Shared types.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>

#include "types.h"

void free_array(struct array_s x) {
    free(x.a);
};

u32 bytes_array(struct array_s *x) {
    return x->length * sizeof(x->a[0]) + sizeof(x->length);
};
