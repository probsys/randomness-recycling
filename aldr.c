/*
  Name:     aldr.c
  Purpose:  Amplified Loaded Dice Roller.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>
#include <string.h>

#include "aldr.h"
#include "uniform.h"

struct aldr_recycle_s preprocess_aldr_recycle(u32* a, u32 n) {
    // assume k <= 31
    u32 m = 0;
    for (u32 i = 0; i < n; ++i) {
        m += a[i];
    }
    u32 k = 32 - __builtin_clz(m) - (0 == (m & (m-1)));
    u32 K = k << 1;
    u64 c = (1ull << K) / m;
    u32 r = (1ull << K) % m;
    u64 *Q = calloc(n, sizeof(u64));

    u32 num_leaves = 0;
    for (u32 i = 0; i < n; ++i) {
        num_leaves += __builtin_popcountll(Q[i] = c * a[i]);
    }

    u32 num_levels = K + 1;
    u32 *breadths = calloc(num_levels, sizeof(u32));
    u32 *leaves_flat = calloc(num_leaves, sizeof(u32));

    u32 location = 0;
    for (u32 j = 0; j <= K; ++j) {
        u64 bit = (1ull << (K - j));
        for (u32 i = 0; i < n; ++i) {
            if (Q[i] & bit) {
                leaves_flat[location] = i;
                ++breadths[j];
                ++location;
            }
        }
    }

    return (struct aldr_recycle_s){
            .length_breadths = num_levels,
            .length_leaves_flat = num_leaves,
            .length_weights = n,
            .reject_weight = r,
            .breadths = breadths,
            .leaves_flat = leaves_flat,
            .weights = Q
        };
}

u32 sample_aldr_recycle(struct aldr_recycle_s* f) {
    u32 num_flips = f->length_breadths - 1;
    while (1) {
        u64 flips = flip_n_from_unif(num_flips);
        if (unlikely(flips >= (1ull << num_flips) - f->reject_weight)) {
            merge_state(flips - (1ull << num_flips) + f->reject_weight, f->reject_weight);
            continue;
        }
        u32 depth = 0;
        u32 location = 0;
        u32 val = 0;
        u32 pos = num_flips;
        for (;;) {
            if (val < f->breadths[depth]) {
                u32 ans = f->leaves_flat[location + val];
                u64 mask = (1ull<<pos) - 1;
                u64 recycle_state = mask & flips;
                u64 recycle_bound = f->weights[ans];
                recycle_state += recycle_bound & mask;
                merge_state(recycle_state, recycle_bound);
                return ans;
            }
            location += f->breadths[depth];
            --pos;
            val = ((val - f->breadths[depth]) << 1) | ((flips >> pos) & 1);
            ++depth;
        }
    }
}

u32 bytes_aldr_recycle(struct aldr_recycle_s *x) {
    return
        sizeof(x->length_breadths)
            + sizeof(x->length_leaves_flat)
            + sizeof(x->length_weights)
            + x->length_breadths * sizeof(x->breadths[0])
            + x->length_leaves_flat * sizeof(x->leaves_flat[0])
            + x->length_weights * sizeof(x->weights[0]);
}

void free_aldr_recycle (struct aldr_recycle_s x) {
    free(x.breadths);
    free(x.leaves_flat);
    free(x.weights);
}


struct fldr_eo_s preprocess_fldr_eo(u32* a, u32 n) {
    // assume k <= 31
    u32 m = 0;
    for (u32 i = 0; i < n; ++i) {
        m += a[i];
    }
    u32 k = 32 - __builtin_clz(m) - (0 == (m & (m-1)));

    u32 num_leaves = 0;
    for (u32 i = 0; i < n; ++i) {
        num_leaves += __builtin_popcount(a[i]);
    }

    u32 num_levels = k + 1;
    u32 *breadths = calloc(num_levels, sizeof(u32));
    u32 *leaves_flat = calloc(num_leaves, sizeof(u32));

    u32 location = 0;
    for(u32 j = 0; j <= k; ++j) {
        u32 bit = (1u << (k - j));
        for (u32 i = 0; i < n; ++i) {
            if (a[i] & bit) {
                leaves_flat[location] = i;
                ++location;
                ++breadths[j];
            }
        }
    }

    u32 *weights = malloc(n * sizeof(u32));
    memcpy(weights, a, n * sizeof(u32));

    return (struct fldr_eo_s){
            .length_breadths = num_levels,
            .length_leaves_flat = num_leaves,
            .length_weights = n,
            .uniform_preprocessed = uniform_preprocess(m),
            .breadths = breadths,
            .leaves_flat = leaves_flat,
            .weights = weights
        };
}

u32 sample_fldr_eo(struct fldr_eo_s* f) {
    u32 num_flips = f->length_breadths - 1;
    u32 depth = 0;
    u32 location = 0;
    u32 val = 0;
    u32 flips = uniform_prediv(&(f->uniform_preprocessed));
    u32 pos = num_flips;
    for (;;) {
        if (val < f->breadths[depth]) {
            u32 ans = f->leaves_flat[location + val];
            u32 mask = (1u<<pos) - 1;
            u32 recycle_state = mask & flips;
            u32 recycle_bound = f->weights[ans];
            recycle_state += recycle_bound & mask;
            // equivalent and maybe faster:
            // recycle_state |= recycle_bound & -(1u<<(pos+1));
            merge_state(recycle_state, recycle_bound);
            return ans;
        }
        location += f->breadths[depth];
        --pos;
        val = ((val - f->breadths[depth]) << 1) | ((flips >> pos) & 1);
        ++depth;
    }
}

u32 bytes_fldr_eo(struct fldr_eo_s *x) {
    return
        sizeof(x->length_breadths)
            + sizeof(x->length_leaves_flat)
            + sizeof(x->length_weights)
            + sizeof(x->uniform_preprocessed)
            + x->length_breadths * sizeof(x->breadths[0])
            + x->length_leaves_flat * sizeof(x->leaves_flat[0])
            + x->length_weights * sizeof(x->weights[0]);
}

void free_fldr_eo(struct fldr_eo_s x) {
    free(x.breadths);
    free(x.leaves_flat);
    free(x.weights);
}
