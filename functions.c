/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include "functions.h"
#include "constants.h"
#include "hypergeom.h"
#include <mpfr.h>


double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j);
double specie_c(double k_perp, double omega);
double specie_h(double k_perp, double omega);


double D(const double k_perp, const double omega)
{
        int p = 0;

        // We start by setting the default precision for MPFR variables based on the value of k_perp. The larger it is the higher the precision required.
        p = 1 + (int) (k_perp / 30);

        mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

        // ToDo: Re-add specie-c
        // return 1 + (specie_c(k_perp, omega) + specie_h(k_perp, omega));

        double r =  1 + (specie_h(k_perp, omega));

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        return r;
}


double specie_c(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, LAMBDA_KAPPA_C_p2, KAPPA_C, OMEGA_CC, RHO_C);
}

double specie_h(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, LAMBDA_KAPPA_H_p2, KAPPA_H, OMEGA_CH, RHO_H);
}


double specie_j(const double k_perp, const double omega, const double lambda_kappa_j_p2, const double kappa_j, const double omega_cj, const double rho_j)
{
        mpfr_t omega_by_omega_cj, two_lambda_j_prime, coeff, h1f2, h2f3, result;
        mpfr_inits(omega_by_omega_cj, two_lambda_j_prime, coeff, h1f2, h2f3, result, (mpfr_ptr) 0);

        calc_omega_by_omega_cj(omega_by_omega_cj, omega, omega_cj);
        calc_two_lambda_j_prime(two_lambda_j_prime, kappa_j, rho_j, k_perp);

        struct coeffs_1f2 c_1f2;
        struct coeffs_2f3 c_2f3;

        init_coeffs(& c_1f2, & c_2f3);

        calc_coeffs_1f2(& c_1f2, kappa_j, omega_by_omega_cj);
        calc_coeffs_2f3(& c_2f3, kappa_j, omega_by_omega_cj);

        mpfr_set_d(result, 1, RND);
        calc_coeff(coeff, omega_by_omega_cj, kappa_j, two_lambda_j_prime);
        hyp1F2(h1f2, c_1f2, two_lambda_j_prime);
        hyp2F3(h2f3, c_2f3, two_lambda_j_prime);

        mpfr_mul(coeff, coeff, h1f2, RND);
        mpfr_add(result, result, coeff, RND);
        mpfr_sub(result, result, h2f3, RND);

        if (FLAG_DENOM)
                mpfr_div_d(result, result, pow(k_perp, 2) * lambda_kappa_j_p2, RND);

        double r = mpfr_get_d(result, RND);
        mpfr_clears(omega_by_omega_cj, two_lambda_j_prime, coeff, h1f2, h2f3, result, (mpfr_ptr) 0);
        clear_coeffs(& c_1f2, & c_2f3);

        return r;
}


/*
 * Define utility functions for calculating the coefficient, pFq coeffs, and two_lambda_j_prime;
 */

void calc_two_lambda_j_prime(mpfr_t result, const double kappa_j, const double rho_j, const double k_perp)
{
        mpfr_set_d(result, k_perp, RND);
        mpfr_mul_d(result, result, rho_j, RND);
        mpfr_pow_ui(result, result, 2, RND);
        mpfr_mul_d(result, result, 2 * (kappa_j - 1.5), RND);
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

        mpfr_set_d(x, kappa_j + 1.5, RND);
        mpfr_add(x, x, omega_by_omega_cj, RND);
        mpfr_gamma(x, x, RND);
        mpfr_div(coeff, coeff, x, RND);

        mpfr_set_d(x, kappa_j + 1.5, RND);
        mpfr_sub(x, x, omega_by_omega_cj, RND);
        mpfr_gamma(x, x, RND);
        mpfr_div(coeff, coeff, x, RND);

        mpfr_set(x, two_lambda_j_prime, RND);
        mpfr_set_d(y, kappa_j + 0.5, RND);
        mpfr_pow(x, x, y, RND);                 // Set x = pow(x,y)
        mpfr_mul(coeff, coeff, x, RND);

        mpfr_clears(pi, csc, x, y, (mpfr_ptr) 0);
        mpfr_free_cache();                              // To clear the creation of the constant pi
}
