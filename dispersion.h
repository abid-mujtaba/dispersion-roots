#include <mpfr.h>


#define NUM_MPFR_VARIABLES 10

double D(double k_perp, double omega);

void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp, mpfr_t * vars);
void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj);
