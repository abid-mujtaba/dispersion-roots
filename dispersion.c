/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include <math.h>
#include <mpfr.h>
#include <gsl/gsl_sf_gamma.h>
#include "dispersion.h"
#include "constants.h"
#include "derived.h"


// Function prototypes
double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j, mpfr_t * vars);
// double specie_j_zero(double kappa_j, double rho_j, struct coeffs_2f3 c_2f3, double lambda_kappa_j_p2);
double specie_c(double k_perp, double omega, mpfr_t * vars);
double specie_h(double k_perp, double omega, mpfr_t * vars);

double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e);
void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp, mpfr_t * vars);
double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_by_n0e);
void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj);

// The following functions are declared here but are defined elsewhere
void calc_first(mpfr_t first, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_second(mpfr_t second, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_third(mpfr_t third, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);


double D(const double k_perp, const double omega)
{
        int p = 0;

        // We start by setting the default precision for MPFR variables based on the value of k_perp. The larger it is the higher the precision required.
        p = 1 + (int) (k_perp / 30);

        mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

        // Now we create dummy variables which will be used throughout the calculations to save time on initialization and clearing at every step of the calculation.
        mpfr_t vars[NUM_MPFR_VARIABLES];

        // Initialize the variables
        for (int i = 0; i < NUM_MPFR_VARIABLES; ++i)
                mpfr_init(vars[i]);


        double r = 1 + (specie_c(k_perp, omega, vars) + specie_h(k_perp, omega, vars));


        // Clear the variables
        for (int i = 0; i < NUM_MPFR_VARIABLES; ++i)
                mpfr_clear(vars[i]);

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        return r;
}


double specie_c(const double k_perp, const double omega, mpfr_t * vars)
{
        return specie_j(k_perp, omega, KAPPA_C, OMEGA_CC, RHO_C, N0C_BY_N0E, vars);
}

double specie_h(const double k_perp, const double omega, mpfr_t * vars)
{
        return specie_j(k_perp, omega, KAPPA_H, OMEGA_CH, RHO_H, N0H_BY_N0E, vars);
}


double specie_j(const double k_perp, const double omega, const double kappa_j, const double omega_cj, const double rho_j, const double n0j_by_n0e, mpfr_t * vars)
{
        // ToDo: Deal with the special case k_perp == 0

        double r;

        mpfr_t result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, term;
        mpfr_inits(result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, term, (mpfr_ptr) 0);

        // Calculate MPFR variables required for the three terms
        mpfr_set_d(kappa, kappa_j, RND);
        calc_omega_by_omega_cj(omega_by_omega_cj, omega, omega_cj);
        calc_two_lambda_j(two_lambda_j, kappa, rho_j, k_perp, vars);
        mpfr_const_pi(pi, RND);

        mpfr_mul(csc, omega_by_omega_cj, pi, RND);    // csc = omega_by_omega_cj * pi
        mpfr_csc(csc, csc, RND);                      // csc = cosec( csc )


        calc_first(result, kappa, omega_by_omega_cj, two_lambda_j, csc, pi, * vars, * (vars + 1), vars + 2);
        calc_second(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi, * vars, * (vars + 1), vars + 2);
        mpfr_sub(result, result, term, RND);         // result = first - second

        calc_third(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi, * vars, * (vars + 1), vars + 2);
        mpfr_sub(result, result, term, RND);         // result -= third


        // Final division
        if (k_perp == 0)
                mpfr_div_d(result, result, lambda_vcj_p2(kappa_j, rho_j, n0j_by_n0e), RND);
        else
                mpfr_div_d(result, result, pow(k_perp, 2) * lambda_vcj_p2(kappa_j, rho_j, n0j_by_n0e), RND);         // result /= pow(k_perp, 2) * lambda_vcj_p2

        r = mpfr_get_d(result, RND);

        mpfr_clears(result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, term, (mpfr_ptr) 0);
        mpfr_free_cache();              // Clear the creation of the constant pi

        return r;
}


void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp, mpfr_t * vars)
{
        mpfr_t * x = vars;              // x equals the first pointer in vars, which has already been initialized

        mpfr_set_d(result, k_perp, RND);                // result = k_perp;
        mpfr_mul_d(result, result, rho_j, RND);         // result *= rho_j;
        mpfr_pow_ui(result, result, 2, RND);            // result = pow(result, 2);

        // We de-reference the pointer x to gain access to the mpfr_t variable it points to
        mpfr_sub_d(*x, kappa_j, 1.5, RND);               // x = kappa - 1.5
        mpfr_mul_ui(*x, *x, 2, RND);                      // x *= 2

        mpfr_mul(result, result, *x, RND);               // result *= x
}


double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return (3 * LAMBDA + 1) * (kappa_j + 1) * gsl_sf_gamma(kappa_j + 1.5) * lambda_kappa_j_p2(kappa_j, rho_j, n0j_by_n0e);
}


double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return (kappa_j - 1.5) * pow(rho_j, 2) / (n0j_by_n0e * (pow(OMEGA_UH_BY_OMEGA_CE, 2) - 1) * (kappa_j - 0.5));
}


void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj)
{
        mpfr_set_d(result, omega, RND);
        mpfr_mul_d(result, result, omega_cj, RND);
}
