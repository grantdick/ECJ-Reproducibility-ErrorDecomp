#ifndef RNG_H
#define RNG_H

#ifdef  __cplusplus
extern "C" {
#endif

    #include <stdio.h>
    #include <stdint.h>

    void seed_rng(const char *const src, uint32_t offset);

    double next_rnd(void);

    double next_rnd_gauss(double mu, double sd);

#ifdef  __cplusplus
}
#endif

#endif
