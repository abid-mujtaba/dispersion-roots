#include <mpfr.h>               // To gain access to the MPFR_RNDN constant, and the mpfr_t type


struct coeffs_1f2 {
        double a1;
        double b1;
        double b2;
};


struct coeffs_2f3 {
        double a1;
        double a2;
        double b1;
        double b2;
        double b3;
};


void hyp1F2(mpfr_t result, const struct coeffs_1f2 c_1f2, const double x);
void hyp2F3(mpfr_t result, const struct coeffs_2f3 c_2f3, const double x);
