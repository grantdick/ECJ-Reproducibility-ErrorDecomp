/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
   numbers. We suggest to use its upper bits for floating-point
   generation, as it is slightly faster than xoshiro256++/xoshiro256**. It
   passes all tests we are aware of except for the lowest three bits,
   which might fail linearity tests (and just those), so if low linear
   complexity is not considered an issue (as it is usually the case) it
   can be used to generate 64-bit outputs, too.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */
#include <stdio.h>
#include <stdlib.h>

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include <sys/time.h>

#include <errno.h>
#include <string.h>

static uint64_t state[4];

void seed_rng(const char *const src, uint32_t offset)
{
    uint8_t i;
    uint64_t x;
    struct timeval t;
    const char *const default_src = "/dev/urandom";
    FILE *rnd_seed;

    if (src == NULL) {
        /* first try random device, and if that cannot be found, then by system time */
        rnd_seed = fopen(default_src, "rb");
        if (rnd_seed == NULL) { /* doesn't look like device exists! */
            gettimeofday(&t, NULL);
            x = t.tv_sec * t.tv_usec;

            /* now initialise the seed with splitmix64 */
            for (i = 0; i < 4; ++i) {
                state[i] = (x += 0x9e3779b97f4a7c15);
                state[i] = (state[i] ^ (state[i] >> 30)) * 0xbf58476d1ce4e5b9;
                state[i] = (state[i] ^ (state[i] >> 27)) * 0x94d049bb133111eb;
                state[i] = state[i] ^ (state[i] >> 31);
            }
        } else {
            if (fread(state, sizeof(uint64_t), 4, rnd_seed) != 4) {
                fprintf(stderr, "Something went wrong with reading seeding from %s\n", default_src);
                exit(EXIT_FAILURE);
            }
            fclose(rnd_seed);
        }
    } else {
        rnd_seed = fopen(src, "rb");

        if (rnd_seed == NULL) {
            fprintf(stderr, "Failed to open %s for random number seeding. Reason: %d (%s)\n", src, errno, strerror(errno));
            exit(EXIT_FAILURE);
        }

        if (fseek(rnd_seed, 4 * offset, SEEK_SET) != 0) {
            fprintf(stderr, "Something went wrong with seeking seed position from %s. Reason: %d (%s)\n", src, errno, strerror(errno));
            exit(EXIT_FAILURE);
        } else if (fread(state, sizeof(uint64_t), 4, rnd_seed) != 4) {
            fprintf(stderr, "Something went wrong with reading seeding from %s\n", src);
            exit(EXIT_FAILURE);
        }

        fclose(rnd_seed);
    }
}

static inline uint64_t rotl(const uint64_t x, uint32_t k)
{
    return (x << k) | (x >> (64 - k));
}

static uint64_t next(uint64_t state[4])
{
    const uint64_t result = state[0] + state[3];

    const uint64_t t = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= t;

    state[3] = rotl(state[3], 45);

    return result;
}

double next_rnd(void)
{
    return (next(state) >> 11) * 0x1.0p-53;
}

double next_rnd_gauss(double mu, double sd)
{
    static bool cached = false;
    static double spare = NAN;

    double u, v, s;

    if (cached) {
        cached = false;
        return mu + sd * spare;
    }

    do {
        u = next_rnd() * 2.0 - 1.0;
        v = next_rnd() * 2.0 - 1.0;
        s = u * u + v * v;
    } while((s >= 1.0) || (s == 0.0));

    s = sqrt(-2.0 * log(s) / s);

    cached = true;
    spare = v * s;

    return mu + sd * u * s;
}
