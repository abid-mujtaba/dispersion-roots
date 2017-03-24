/*
 * Define the prototypes of functions defined in functions.c
 */

double D(double k_perp, double omega);

double calculate_two_lambda_j_prime(const double kappa_j, const double rho_j, const double k_perp);
double coefficient(const double omega_by_omega_cj, const double kappa_j, const double two_lambda_j_prime);

#define FLAG_DENOM 1
