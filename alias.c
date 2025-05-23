/*
  Name:     alias.c
  Purpose:  Exact weighted alias method, translated from Rust to C.
  Author:   CMU Probabilistic Computing Systems Lab
  Copyright (C) 2025 CMU Probabilistic Computing Systems Lab, All Rights Reserved.

  Released under Apache 2.0; refer to LICENSE.txt
*/

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

#include "uniform.h"
#include "types.h"
#include "alias.h"

void free_weighted_alias(struct weighted_alias_s x) {
    free(x.aliases);
    free(x.no_alias_odds);
}

struct Aliases {
    u32 *aliases;
    u32 smalls_head;
    u32 bigs_head;
};

/// This struct is designed to contain three data structures at once,
/// sharing the same memory. More precisely it contains two linked lists
/// and an alias map, which will be the output of this method. To keep
/// the three data structures from getting in each other's way, it must
/// be ensured that a single index is only ever in one of them at the
/// same time.
struct Aliases aliases_new(u32 n) {
    return (struct Aliases) {
        .aliases = calloc(n, sizeof(u32)),
        .smalls_head = UINT32_MAX,
        .bigs_head = UINT32_MAX
    };
}

void push_small(struct Aliases *aliases, u32 idx) {
    aliases->aliases[idx] = aliases->smalls_head;
    aliases->smalls_head = idx;
}

void push_big(struct Aliases *aliases, u32 idx) {
    aliases->aliases[idx] = aliases->bigs_head;
    aliases->bigs_head = idx;
}

u32 pop_small(struct Aliases *aliases) {
    u32 idx = aliases->smalls_head;
    aliases->smalls_head = aliases->aliases[idx];
    return idx;
}

u32 pop_big(struct Aliases *aliases) {
    u32 idx = aliases->bigs_head;
    aliases->bigs_head = aliases->aliases[idx];
    return idx;
}

bool smalls_is_empty(struct Aliases *aliases) {
    return aliases->smalls_head == UINT32_MAX;
}

bool bigs_is_empty(struct Aliases *aliases) {
    return aliases->bigs_head == UINT32_MAX;
}

void set_alias(struct Aliases *aliases, u32 idx, u32 alias) {
    aliases->aliases[idx] = alias;
}

/// Creates a new [`WeightedAliasIndex`].
///
/// Returns an error if:
/// - The vector is empty.
/// - The vector is longer than `u32::MAX`.
/// - For any weight `w`: `w < 0` or `w > max` where `max = W::MAX /
///   weights.len()`.
/// - The sum of weights is zero.
struct weighted_alias_s preprocess_weighted_alias(int* a, int n) {
    assert(n > 0);
    assert(n < UINT32_MAX);
    u32 max_weight_size = UINT32_MAX / n;
    for (u32 i = 0; i < n; ++i) {
        assert(0 <= a[i]);
        assert(a[i] <= max_weight_size);
    }

    // The sum of weights will represent 100% of no alias odds.
    u32 weight_sum = 0;
    for (u32 i = 0; i < n; ++i) {
        weight_sum += a[i];
    }
    assert(weight_sum >= 0);

    u32 *no_alias_odds = calloc(n, sizeof(u32));
    for (u32 i = 0; i < n; ++i) {
        no_alias_odds[i] = a[i] * n;
    }

    struct Aliases aliases = aliases_new(n);

    // Split indices into those with small weights and those with big weights.
    for (u32 i = 0; i < n; ++i) {
        if (no_alias_odds[i] < weight_sum) {
            push_small(&aliases, i);
        } else {
            push_big(&aliases, i);
        }
    }

    // Build the alias map by finding an alias with big weight for each index with
    // small weight.
    while (!smalls_is_empty(&aliases) && !bigs_is_empty(&aliases)) {
        u32 small = pop_small(&aliases);
        u32 big = pop_big(&aliases);
        set_alias(&aliases, small, big);
        no_alias_odds[big] -= weight_sum - no_alias_odds[small];
        if (no_alias_odds[big] < weight_sum) {
            push_small(&aliases, big);
        } else {
            push_big(&aliases, big);
        }
    }

    // The remaining indices should have no alias odds of about 100%. This is due to
    // numeric accuracy. Otherwise they would be exactly 100%.
    while (!smalls_is_empty(&aliases)) {
        no_alias_odds[pop_small(&aliases)] = weight_sum;
    }
    while (!bigs_is_empty(&aliases)) {
        no_alias_odds[pop_big(&aliases)] = weight_sum;
    }

    return (struct weighted_alias_s) {
        .length = n,
        .weight_sum = weight_sum,
        .aliases = aliases.aliases,
        .no_alias_odds = no_alias_odds
    };
}

int bytes_weighted_alias(struct weighted_alias_s *x) {
    return
        x->length * sizeof(x->aliases[0])
            + x->length * sizeof(x->no_alias_odds[0])
            + sizeof(x->length)
            + sizeof(x->weight_sum);
}

u32 sample_weighted_alias_recycle(struct weighted_alias_s *x) {
    u32 uniform_index = uniform_eo(x->length);
    if (bernoulli_eo(x->no_alias_odds[uniform_index], x->weight_sum)) {
        return uniform_index;
    } else {
        return x->aliases[uniform_index];
    }
}

struct weighted_alias_eo_s preprocess_weighted_alias_eo(int* a, int n) {
    struct weighted_alias_s wai = preprocess_weighted_alias(a, n);

    u64 *cumulative_sums = calloc(wai.length, sizeof(u64));
    for (u32 i = 0; i < wai.length; ++i) {
        cumulative_sums[i] = wai.no_alias_odds[i];
    }
    u64 *offsets = calloc(wai.length, sizeof(u64));
    for (u32 i = 0; i < wai.length; ++i) {
        if (wai.aliases[i] != UINT32_MAX) {
            // might underflow but doesn't matter:
            offsets[i] = cumulative_sums[wai.aliases[i]] - wai.no_alias_odds[i];
            cumulative_sums[wai.aliases[i]] += wai.weight_sum - wai.no_alias_odds[i];
        }
    }
    free(cumulative_sums);
    u32 *weights = malloc(wai.length * sizeof(u32));
    memcpy(weights, a, wai.length * sizeof(u32));

    return (struct weighted_alias_eo_s) {
        .length = wai.length,
        .weight_sum = wai.weight_sum,
        .weights = weights,
        .aliases = wai.aliases,
        .no_alias_odds = wai.no_alias_odds,
        .offsets = offsets
    };
}

int bytes_weighted_alias_eo(struct weighted_alias_eo_s *x) {
    return
        x->length * sizeof(x->aliases[0])
            + x->length * sizeof(x->no_alias_odds[0])
            + x->length * sizeof(x->weights[0])
            + x->length * sizeof(x->offsets[0])
            + sizeof(x->length)
            + sizeof(x->weight_sum);
}

u32 sample_weighted_alias_eo(struct weighted_alias_eo_s *x) {
    u64 uniform_index = uniform_eo((u64)x->length * (u64)x->weight_sum);
    u64 uniform_weight = uniform_index / x->length;
    uniform_index %= x->length;
    u64 no_alias_odds = x->no_alias_odds[uniform_index];
    if (uniform_weight < no_alias_odds) {
        merge_state(uniform_weight, (u64)x->weights[uniform_index] * (u64)x->length);
        return uniform_index;
    } else {
        merge_state(uniform_weight + x->offsets[uniform_index], (u64)x->weights[x->aliases[uniform_index]] * (u64)x->length);
        return x->aliases[uniform_index];
    }
}

void free_weighted_alias_eo(struct weighted_alias_eo_s x) {
    free(x.weights);
    free(x.aliases);
    free(x.no_alias_odds);
    free(x.offsets);
}
