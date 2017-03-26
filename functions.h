/*
 * Define the prototypes of functions defined in functions.c
 */

 #include <mpfr.h>


double D(double k_perp, double omega);

double calc_two_lambda_j_prime(const double kappa_j, const double rho_j, const double k_perp);
struct coeffs_1f2 calc_coeffs_1f2(const double kappa_j, const mpfr_t omega_by_omega_cj);
struct coeffs_2f3 calc_coeffs_2f3(const double kappa_j, const mpfr_t omega_by_omega_cj);
void calc_coeff(mpfr_t result, const mpfr_t omega_by_omega_cj, const double kappa_j, const double two_lambda_j_prime);

#define FLAG_DENOM 1
