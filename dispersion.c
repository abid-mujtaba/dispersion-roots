/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include <mpfr.h>
#include <gsl/gsl_sf_gamma.h>
#include "dispersion.h"
#include "henning.h"            // For lambda_kappa_j_p2(,,)
#include "constants.h"
#include "hypergeom.h"


// Function prototypes
double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j);
double specie_j_zero(double kappa_j, double rho_j, struct coeffs_2f3 c_2f3, double lambda_kappa_j_p2);
double specie_c(double k_perp, double omega);
double specie_h(double k_perp, double omega);

double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e);
void calc_two_lambda_j_prime(mpfr_t result, const double kappa_j, const double rho_j, const double k_perp);

void calc_coeffs_1f2(struct coeffs_1f2 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj);
void calc_coeffs_2f3(struct coeffs_2f3 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj);

void calc_coeff(mpfr_t result, const mpfr_t omega_by_omega_cj, const double kappa_j, const mpfr_t two_lambda_j_prime);
void calc_unnorm_coeff(mpfr_t coeff, const mpfr_t omega_by_omega_cj, const double kappa_j, const mpfr_t two_lambda_j_prime);


double D(const double k_perp, const double omega)
{
        // Testing using printf statements
        printf("\nlambda_vcj_p2 = %.17g", lambda_vcj_p2(KAPPA_H, RHO_H, N0H_BY_N0E));


        int p = 0;

        // We start by setting the default precision for MPFR variables based on the value of k_perp. The larger it is the higher the precision required.
        p = 1 + (int) (k_perp / 30);

        mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

        double r = 1 + (specie_c(k_perp, omega) + specie_h(k_perp, omega));

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        return r;
}


double specie_c(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, KAPPA_C, OMEGA_CC, RHO_C, N0C_BY_N0E);
}

double specie_h(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, KAPPA_H, OMEGA_CH, RHO_H, N0H_BY_N0E);
}


double specie_j(const double k_perp, const double omega, const double kappa_j, const double omega_cj, const double rho_j, const double n0j_by_n0e)
{
        // ToDo: Deal with the special case k_perp == 0

        return 0;
}

/*
 * Calculate specie_j for the special case of k_perp = 0.
 * In this case the only term that survives is the first term of 2F3, all other terms going to zero.
 */
double specie_j_zero(double kappa_j, double rho_j, struct coeffs_2f3 c, double lambda_kappa_j_p2)
{
        mpfr_t q;
        mpfr_init(q);

        mpfr_set_d(q, 1, RND);
        mpfr_mul(q, q, c.a1, RND);
        mpfr_mul(q, q, c.a2, RND);
        mpfr_div(q, q, c.b1, RND);
        mpfr_div(q, q, c.b2, RND);
        mpfr_div(q, q, c.b3, RND);

        // Multiply the coefficient of k_perp^2 in two_lambda_j_prime
        mpfr_mul_d(q, q, 2 * (kappa_j - 1.5), RND);
        mpfr_mul_d(q, q, pow(rho_j, 2), RND);

        // Divide by the lambda_kappa_j_p2 in the denominator with k_perp^2 (final part of calculation of specie_j)
        mpfr_div_d(q, q, lambda_kappa_j_p2, RND);

        double r = mpfr_get_d(q, RND);
        mpfr_clear(q);

        return -r;
}


/*
 * Define utility functions for calculating the coefficient, pFq coeffs, and two_lambda_j_prime;
 */

void calc_two_lambda_j_prime(mpfr_t result, const double kappa_j, const double rho_j, const double k_perp)
{
        mpfr_set_d(result, k_perp, RND);                // result = k_perp;
        mpfr_mul_d(result, result, rho_j, RND);         // result *= rho_j;
        mpfr_pow_ui(result, result, 2, RND);            // result = pow(result, 2);
        mpfr_mul_d(result, result, 2 * (kappa_j - 1.5), RND);  // result *= 2 * (kappa_j - 1.5);
}


void calc_coeffs_1f2(struct coeffs_1f2 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj)
{
        mpfr_set_d(c->a1, kappa_j + 1, RND);
        mpfr_set_d(c->b1, kappa_j + 1.5, RND);
        mpfr_set_d(c->b2, kappa_j + 1.5, RND);

        mpfr_add(c->b1, c->b1, omega_by_omega_cj, RND);
        mpfr_sub(c->b2, c->b2, omega_by_omega_cj, RND);
}


void calc_coeffs_2f3(struct coeffs_2f3 * const c, const double kappa_j, const mpfr_t omega_by_omega_cj)
{
        mpfr_set_d(c->a1, 1, RND);
        mpfr_set_d(c->a2, 0.5, RND);
        mpfr_set_d(c->b1, 0.5 - kappa_j, RND);
        mpfr_set_d(c->b2, 1, RND);
        mpfr_set_d(c->b3, 1, RND);

        mpfr_add(c->b2, c->b2, omega_by_omega_cj, RND);
        mpfr_sub(c->b3, c->b3, omega_by_omega_cj, RND);
}


void calc_coeff(mpfr_t coeff, const mpfr_t omega_by_omega_cj, const double kappa_j, const mpfr_t two_lambda_j_prime)
{
        mpfr_t x;
        mpfr_init(x);

        calc_unnorm_coeff(coeff, omega_by_omega_cj, kappa_j, two_lambda_j_prime);

        mpfr_set_d(x, kappa_j + 1.5, RND);
        mpfr_sub(x, x, omega_by_omega_cj, RND);
        mpfr_gamma(x, x, RND);
        mpfr_div(coeff, coeff, x, RND);

        mpfr_clear(x);
}


/* Coeff calculating without dividing by the gamma function of b2.
 * This is used when one uses norm_1F2 which contains the above-mentioned division within itself.
 */
void calc_unnorm_coeff(mpfr_t coeff, const mpfr_t omega_by_omega_cj, const double kappa_j, const mpfr_t two_lambda_j_prime)
{
        mpfr_t pi, csc, x, y;
        mpfr_inits(pi, csc, x, y, (mpfr_ptr) 0);

        mpfr_const_pi(pi, RND);              // Calculate pi and store it in the variable
        mpfr_sqrt(coeff, pi, RND);         // Calculate the square root of the second argument and store it in the first

        mpfr_mul(coeff, coeff, omega_by_omega_cj, RND);

        mpfr_set_d(x, kappa_j + 1, RND);                // Store 'kappa_j + 1' in x as mpfr value
        mpfr_gamma(x, x, RND);                          // Calculate gamma(x) and store in x
        mpfr_mul(coeff, coeff, x, RND);                 // Multiply coeff with this value (gamma(x))

        mpfr_set_d(x, 0.5 - kappa_j, RND);
        mpfr_gamma(x, x, RND);
        mpfr_mul(coeff, coeff, x, RND);

        mpfr_set(x, omega_by_omega_cj, RND);
        mpfr_mul(x, x, pi, RND);
        mpfr_csc(csc, x, RND);
        mpfr_mul(coeff, coeff, csc, RND);

        mpfr_set(x, two_lambda_j_prime, RND);
        mpfr_set_d(y, kappa_j + 0.5, RND);
        mpfr_pow(x, x, y, RND);                 // Set x = pow(x,y)
        mpfr_mul(coeff, coeff, x, RND);

        mpfr_set_d(x, kappa_j + 1.5, RND);              // Divide by Gamma(c.b1) but NOT c.b2
        mpfr_add(x, x, omega_by_omega_cj, RND);
        mpfr_gamma(x, x, RND);
        mpfr_div(coeff, coeff, x, RND);

        mpfr_clears(pi, csc, x, y, (mpfr_ptr) 0);
        mpfr_free_cache();                              // To clear the creation of the constant pi
}


double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return (3 * LAMBDA + 1) * (kappa_j + 1) * gsl_sf_gamma(kappa_j + 1.5) * lambda_kappa_j_p2(kappa_j, rho_j, n0j_by_n0e);
}
