/*
 * Define the prototypes of functions defined in functions.c
 */

 #include <mpfr.h>
 #include "hypergeom.h"


double D(double k_perp, double omega);

void calc_two_lambda_j_prime(mpfr_t result, const double kappa_j, const double rho_j, const double k_perp);
void calc_coeffs_1f2(struct coeffs_1f2 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj);
void calc_coeffs_2f3(struct coeffs_2f3 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj);
void calc_coeff(mpfr_t result, const mpfr_t omega_by_omega_cj, const double kappa_j, const mpfr_t two_lambda_j_prime);

#define FLAG_DENOM 1
