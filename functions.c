/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include "functions.h"
#include "constants.h"
#include "hypergeom.h"
#include <gsl/gsl_math.h>       // Defines M_SQRTPI
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_gamma.h>
#include <mpfr.h>


double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j);
double specie_c(double k_perp, double omega);
double specie_h(double k_perp, double omega);


double D(const double k_perp, const double omega)
{
        // We start by setting the default precision for MPFR variables
        mpfr_set_default_prec(PRECISION);

        // ToDo: Re-add specie-c
        // return 1 + (specie_c(k_perp, omega) + specie_h(k_perp, omega));

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        return 1 + (specie_h(k_perp, omega));
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
        const double two_lambda_j_prime = calc_two_lambda_j_prime(kappa_j, rho_j, k_perp);
        const double omega_by_omega_cj = omega / omega_cj;

        double coeff = calc_coeff(omega_by_omega_cj, kappa_j, two_lambda_j_prime);

        struct coeffs_1f2 c_1f2 = calc_coeffs_1f2(kappa_j, omega_by_omega_cj);
        struct coeffs_2f3 c_2f3 = calc_coeffs_2f3(kappa_j, omega_by_omega_cj);

        double result = 1;
        result += coeff * hyp1F2(c_1f2, two_lambda_j_prime);
        result -= hyp2F3(c_2f3, two_lambda_j_prime);

        // printf("\ncoeff = %.17g", coeff);
        // printf("\n1F2 = %.17g", hyp1F2(c_1f2, two_lambda_j_prime));
        // printf("\n3F2 = %.17g", hyp2F3(c_2f3, two_lambda_j_prime));
        // printf("\nresult (no denom) = %.17g", result);

        if (FLAG_DENOM)
                result /= (pow(k_perp, 2) * lambda_kappa_j_p2);

        return result;
}


/*
 * Define utility functions for calculating the coefficient, pFq coeffs, and two_lambda_j_prime;
 */

double calc_two_lambda_j_prime(const double kappa_j, const double rho_j, const double k_perp)
{
        return 2 * (kappa_j - 1.5) * pow(k_perp * rho_j, 2);
}


struct coeffs_1f2 calc_coeffs_1f2(const double kappa_j, const double omega_by_omega_cj)
{
        struct coeffs_1f2 c;

        c.a1 = kappa_j + 1;
        c.b1 = kappa_j + 1.5 + omega_by_omega_cj;
        c.b2 = kappa_j + 1.5 - omega_by_omega_cj;

        return c;
}


struct coeffs_2f3 calc_coeffs_2f3(const double kappa_j, const double omega_by_omega_cj)
{
        struct coeffs_2f3 c;

        c.a1 = 1;
        c.a2 = 0.5;
        c.b1 = 0.5 - kappa_j;
        c.b2 = 1 + omega_by_omega_cj;
        c.b3 = 1 - omega_by_omega_cj;

        return c;
}


double calc_coeff(const double omega_by_omega_cj, const double kappa_j, const double two_lambda_j_prime)
{
        double r;

        mpfr_t coeff, pi, csc;
        mpfr_inits(coeff, pi, csc, (mpfr_ptr) 0);

        mpfr_const_pi(pi, MPFR_RNDN);              // Calculate pi and store it in the variable
        mpfr_sqrt(coeff, pi, MPFR_RNDN);         // Calculate the square root of the second argument and store it in the first

        mpfr_mul_d(coeff, coeff, omega_by_omega_cj, MPFR_RNDN);

        // NOTE: We have NOT used MPFR for calculating these gamma function values. Loss of precision.
        mpfr_mul_d(coeff, coeff, gsl_sf_gamma(kappa_j + 1) * gsl_sf_gamma(0.5 - kappa_j), MPFR_RNDN);

        mpfr_set_d(csc, omega_by_omega_cj, MPFR_RNDN);
        mpfr_mul(csc, csc, pi, MPFR_RNDN);
        mpfr_csc(csc, csc, MPFR_RNDN);

        mpfr_mul(coeff, coeff, csc, MPFR_RNDN);
        mpfr_div_d(coeff, coeff, gsl_sf_gamma(kappa_j + 1.5 + omega_by_omega_cj) * gsl_sf_gamma(kappa_j + 1.5 - omega_by_omega_cj), RND);
        mpfr_mul_d(coeff, coeff, pow(two_lambda_j_prime, kappa_j + 0.5), RND);

        r = mpfr_get_d(coeff, RND);

        mpfr_clears(coeff, pi, (mpfr_ptr) 0);
        mpfr_free_cache();                              // To clear the creation of the constant pi

        return r;
}
