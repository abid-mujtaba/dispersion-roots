/*
 * Define the constants used externally.
 */

#include "constants.h"
#include <math.h>
#include <mpfr.h>


// This function provides a centralized mechanism for calculating lambda_kappa_j_p2 from the relevant specie variables
double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return (kappa_j - 1.5) * pow(rho_j, 2) / (n0j_by_n0e * (pow(OMEGA_UH_BY_OMEGA_CE, 2) - 1) * (kappa_j - 0.5));
}


void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj)
{
        mpfr_set_d(result, omega, RND);
        mpfr_mul_d(result, result, omega_cj, RND);
}
