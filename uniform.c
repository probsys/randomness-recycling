/*
  Name:     uniform.c
  Purpose:  Generating uniform pseudo-random numbers.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>
#include <sys/random.h>

#include "uniform.h"

const u32 flip_k = 64;
u64 flip_word = 0;
u32 flip_pos = 0;

void refill(void) {
    getrandom(&flip_word, sizeof(flip_word), 0);
    flip_pos = flip_k;
}

void check_refill(void) {
    if (flip_pos == 0) {
        refill();
    }
}

u64 flip_n(u32 n) {
    check_refill();
    u32 num_bits_extract = min(n, flip_pos);
    flip_pos -= num_bits_extract;
    u64 b = (flip_word >> flip_pos) & (UINT64_MAX >> (64 - num_bits_extract));
    if (num_bits_extract != n) {
        refill();
        num_bits_extract = n - num_bits_extract;
        b <<= num_bits_extract;
        flip_pos -= num_bits_extract;
        b |= (flip_word >> flip_pos) & (UINT64_MAX >> (64 - num_bits_extract));
    }
    return b;
}

// unif_state ~ unif[0, unif_bound)
u64 unif_state = 0;
u64 unif_bound = 1;

void check_refill_uniform() {
    // Update unif_state and unif_bound so that
    // unif_bound >= (1<<63),
    // while retaining
    // unif_state ~ unif[0, unif_bound).
    u32 num_bits_extract = __builtin_clzll(unif_bound);
    if (num_bits_extract >= 8) {
        unif_bound <<= num_bits_extract;
        unif_state <<= num_bits_extract;
        unif_state |= flip_n(num_bits_extract);
    }
}

void merge_state(u64 state, u64 bound) {
    // Input state and bound must be
    // independent of unif_state and unif_bound and satisfy
    // state ~ unif[0, bound).
    // Merge them into unif_state and unif_bound,
    // retaining unif_state ~ unif[0, unif_bound).
    unif_bound *= bound;
    unif_state = unif_state * bound + state;
}

void merge_state_bits(u64 state, u64 n) {
    // Specialize merge_state for n-bit states.
    unif_bound <<= n;
    unif_state = (unif_state << n) | state;
}

u64 uniform_eo(u64 n) {
    // Input positive integer n should be (much) smaller than 1<<63.
    // Output is distributed as unif[0, n),
    // while unif_state is independent of the output and retains
    // unif_state ~ unif[0, unif_bound).
    check_refill_uniform();
    u64 q_state = unif_state / n;
    u64 r_state = unif_state % n;
    u64 q_bound = unif_bound / n;
    u64 r_bound = unif_bound % n;
    // Discard information of bernoulli(r_bound, unif_bound)
    // to split into two branches.
    if (likely(q_state < q_bound)) {
        // q_state ~ unif[0, q_bound)
        // r_state ~ unif[0, n)
        // q_state and r_state are independent
        unif_state = q_state;
        unif_bound = q_bound;
        return r_state;
    } else {
        // q_state = q_bound
        // r_state ~ unif[0, r_bound)
        unif_state = r_state;
        unif_bound = r_bound;
        return uniform_eo(n);
    }
}

u64 flip_n_from_unif(u32 n) {
    // Specialize uniform_eo to use bit shifts, not division,
    // for n uniform bits.
    // Use this instead of the random bit source directly
    // if you plan to recycle randomness, to avoid overflow.
    check_refill_uniform();
    u64 q_state = unif_state >> n;
    u64 r_state = unif_state & ((1ull << n) - 1);
    u64 q_bound = unif_bound >> n;
    u64 r_bound = unif_bound & ((1ull << n) - 1);
    if (likely(q_state < q_bound)) {
        unif_state = q_state;
        unif_bound = q_bound;
        return r_state;
    } else {
        unif_state = r_state;
        unif_bound = r_bound;
        return flip_n_from_unif(n);
    }
}

u32 uniform_u32_from_unif() {
    // Specialize uniform_eo to use bit shifts, not division,
    // for the case of n = 1<<32.
    // Use this instead of the random bit source directly
    // if you plan to recycle randomness, to avoid overflow.
    check_refill_uniform();
    u32 q_state = unif_state >> 32;
    u32 r_state = unif_state;
    u32 q_bound = unif_bound >> 32;
    u32 r_bound = unif_bound;
    if (likely(q_state < q_bound)) {
        unif_state = q_state;
        unif_bound = q_bound;
        return r_state;
    } else {
        unif_state = r_state;
        unif_bound = r_bound;
        return uniform_u32_from_unif();
    }
}

struct uniform_preprocessed_s uniform_preprocess(u32 m) {
    // 1 < m < 2^32
    u64 numerator = 1ull << 32;
    u32 quotient = numerator / m;
    u32 remainder = numerator % m;
    u32 not_remainder = ~ remainder;
    u64 inverse = __UINT64_MAX__ / m;
    u64 inverse_remainder = __UINT64_MAX__ % m;
    if (inverse_remainder == m - 1) ++inverse;
    return (struct uniform_preprocessed_s) {
        .num_outcomes = m,
        .quotient = quotient,
        .not_remainder = not_remainder,
        .inverse = inverse
    };
}

u32 uniform_prediv(struct uniform_preprocessed_s *x) {
    // Compute and recycle uniform, with precomputed divisions.
    u32 u = uniform_u32_from_unif();
    u64 unifm_rem = ((u64) u) * x->num_outcomes;
    u32 unifm = unifm_rem >> 32;
    u32 rem = unifm_rem;
    if (unlikely(rem > x->not_remainder)) {
        // Don't bother trying to recycle the remainder
        return uniform_prediv(x);
    } else {
        // Compute ceiling of (1<<32) * (unifm / m), so
        // u-lower_bound ~ unif[0, x->quotient)
        // unifm ~ unif[0, m)
        // u-lower_bound and unifm are independent
        u32 lower_bound = (x->inverse * unifm) >> 32;
        merge_state(u - lower_bound, x->quotient);
        return unifm;
    }
}

bool bernoulli_eo_2div(u32 numer, u32 denom) {
    u32 unif = uniform_eo(denom);
    if (unif < numer) {
        merge_state(unif, numer);
        return 1;
    } else {
        merge_state(unif - numer, denom - numer);
        return 0;
    }
}

bool bernoulli_eo(u32 numer, u32 denom) {
    check_refill_uniform();
    u64 q_bound = unif_bound / denom;
    u64 r_bound = unif_bound % denom;
    u64 true_bound = q_bound * numer;
    if (unif_state < true_bound) {
        unif_bound = true_bound;
        return 1;
    }
    u64 full_bound = q_bound * denom;
    if (likely(unif_state < full_bound)) {
        unif_state -= true_bound;
        unif_bound = full_bound - true_bound;
        return 0;
    }
    unif_state -= full_bound;
    unif_bound = r_bound;
    return bernoulli_eo(numer, denom);
}
