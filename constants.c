/*
 * Define the constants used externally.
 */

#include "constants.h"
#include <mpfr.h>


void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj)
{
        mpfr_set_d(result, omega, RND);
        mpfr_mul_d(result, result, omega_cj, RND);
}
