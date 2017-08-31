#ifndef HYPERGEOM_H
#define HYPERGEOM_H

#include <mpfr.h>               // To gain access to the MPFR_RNDN constant, and the mpfr_t type


struct coeffs_1f2 {
        mpfr_t * a1;
        mpfr_t * b1;
        mpfr_t * b2;
};


struct coeffs_2f2 {
        mpfr_t * a1;
        mpfr_t * a2;
        mpfr_t * b1;
        mpfr_t * b2;
};


struct coeffs_2f3 {
        mpfr_t *  a1;
        mpfr_t *  a2;
        mpfr_t *  b1;
        mpfr_t *  b2;
        mpfr_t *  b3;
};


void hyp2F2(mpfr_t result, const struct coeffs_2f2 c_2f2, const mpfr_t x);
void hyp2F3(mpfr_t result, const struct coeffs_2f3 c_2f3, const mpfr_t x);

// Calculate the hypergeom functions divided by the gamma function for the last
// b_i coefficient
void norm_hyp2F3(mpfr_t result, const struct coeffs_2f3 c_2f3, const mpfr_t x);

#endif
