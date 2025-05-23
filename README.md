# Randomness Recycling

This repository contains a reference implementation in C of
random sampling algorithms for discrete probability distributions.
These samplers use a technique called randomness recycling to
reduce the expected amortized entropy consumption over many samples.

## Installation

The library can be built by running

```sh
make all
```

This command creates the following artifacts in the `build/` directory:

| Path                  | Description                                                   |
| --------------------- | ------------------------------------------------------------- |
| `build/bin/sample_rr` | Executable for command line interface to randomness recycling |
| `build/include`       | Header files for C programs that use randomness recycling     |
| `build/lib/librr.a`   | Static library for C programs that use randomness recycling   |

## Usage

The following code from [examples/example.c](examples/example.c)
shows how to use randomness recycling to repeatedly sample from a distribution
defined by a list of nonnegative integer weights, using either
inversion sampling via CDF binary search or lookup table,
the alias method,
the Fast Loaded Dice Roller, or
the Amplified Loaded Dice Roller.

```c
#include <stdlib.h>
#include <stdio.h>

#include "aldr.h"
#include "alias.h"
#include "binarysearch.h"
#include "lookup.h"
#include "uniform.h"
#include "types.h"

int main(int argc, char **argv) {
    uint32_t num_samples_each = 18;
    uint32_t num_samplers = 5;
    uint32_t num_samples = num_samples_each * num_samplers;
    uint32_t *samples = calloc(num_samples, sizeof(*samples));

    // Generate a uniform sample with randomness recycling.
    uint64_t U = uniform_eo(50);
    printf("uniform sample: %ld\n", U);

    // Generate nonuniform samples with randomness recycling.
    uint32_t distribution[5] = { 1, 1, 2, 3, 2 };
    struct array_s s_cdf = preprocess_cdf(distribution, 5);
    struct lookup_eo_s s_lookup = preprocess_lookup_eo(distribution, 5);
    struct weighted_alias_eo_s s_alias = preprocess_weighted_alias_eo(distribution, 5);
    struct fldr_eo_s s_fldr = preprocess_fldr_eo(distribution, 5);
    struct aldr_recycle_s s_aldr = preprocess_aldr_recycle(distribution, 5);
    for (uint32_t i = 0; i < num_samples_each; ++i) {
        samples[i] = sample_cdf_eo(&s_cdf);
        samples[i+num_samples_each] = sample_lookup_eo(&s_lookup);
        samples[i+2*num_samples_each] = sample_weighted_alias_eo(&s_alias);
        samples[i+3*num_samples_each] = sample_fldr_eo(&s_fldr);
        samples[i+4*num_samples_each] = sample_aldr_recycle(&s_aldr);
    }

    printf("nonuniform samples: ");
    for (u32 i = 0; i < num_samples; ++i) {
        printf("%d ", samples[i]);
    }
    printf("\n");

    free(samples);
    free_array(s_cdf);
    free_lookup_eo(s_lookup);
    free_weighted_alias_eo(s_alias);
    free_fldr_eo(s_fldr);
    free_aldr_recycle(s_aldr);

    return 0;
}
```

## Usage (Command Line Interface)

The executable in `build/bin/sample_rr` has the following command line interface:

```
usage: ./build/bin/sample_rr <sampler (cdf|lookup|alias|fldr|aldr)> <num_samples> <distribution>
```

where `<num_samples>` is an integer denoting the number of samples to draw,
satisfying `0 <= num_samples <= 2147483647`;
and `<distribution>` is a space-separated list of positive integer weights
for the desired discrete distribution,
with the total number of elements bounded as `0 < n <= 2147483647`,
and the sum bounded as `0 < m <= 2147483647`.

For example, to generate 90 samples from { 1, 1, 2, 3, 2 } using the alias method,
run the following:

```sh
./build/bin/sample_rr alias 90 1 1 2 3 2
```

To generate 9000 samples from { 1, 1, 2, 3, 2 } using a lookup table
and count them as a histogram, run the following:

```sh
./build/bin/sample_rr lookup 9000 1 1 2 3 2 | tr -d '\n' | tr ' ' '\n' | sort | uniq -c
```
