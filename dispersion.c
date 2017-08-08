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
#include "infinite_kappa.h"


// Function prototypes
double specie(double k_perp, double omega, struct Constants * const cj, mpfr_t * vars);

// The following functions are declared here but are defined elsewhere
void calc_first(mpfr_t first, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_second(mpfr_t second, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_third(mpfr_t third, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);


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


        // Fetch the constants for each specie
        struct Constants cc, ch;
        get_constants_c(& cc);
        get_constants_h(& ch);


        // Calculate the result by adding the contributions from both the hot and cold species
        // double r = 1 + (specie(k_perp, omega, & cc, vars) + specie(k_perp, omega, & ch, vars));
        double r = 1 + specie_kappa_infinity(k_perp, omega, & ch, vars);


        // Clear the variables
        for (int i = 0; i < NUM_MPFR_VARIABLES; ++i)
                mpfr_clear(vars[i]);

        // Clear variables inside the constants
        clear_constants(& cc);
        clear_constants(& ch);

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        if (isnan(r))
            fprintf(stderr, "\nWarning - D(%.2f, %.2f) = NAN", k_perp, omega);

        return r;
}


double specie(const double k_perp, const double omega, struct Constants * const c, mpfr_t * vars)
{
        double r;

        mpfr_t result, term;
        mpfr_inits(result, term, (mpfr_ptr) 0);

        // Calculate MPFR variables required for the three terms
        calc_omega_by_omega_cj(c->omega_by_omega_c, omega, c->omega_c);
        calc_two_lambda_j(c->two_lambda, c->kappa, c->rho, k_perp, vars);

        mpfr_mul(c->csc, c->omega_by_omega_c, c->pi, RND);    // csc = omega_by_omega_cj * pi
        mpfr_csc(c->csc, c->csc, RND);                        // csc = cosec( csc )


        calc_first(result, c, * vars, * (vars + 1), vars + 2);
        calc_second(term, c, * vars, * (vars + 1), vars + 2);
        mpfr_sub(result, result, term, RND);         // result = first - second

        calc_third(term, c, * vars, * (vars + 1), vars + 2);
        mpfr_sub(result, result, term, RND);         // result -= third


        // Final division
        mpfr_div(result, result, c->lambda_vc_p2, RND);

        if (k_perp != 0)            // If k_perp = 0 then the result has already been calculated by taking the limit and cancelling out the k_perp^2 term in the denominator
                mpfr_div_d(result, result, pow(k_perp, 2), RND);

        r = mpfr_get_d(result, RND);

        mpfr_clears(result, term, (mpfr_ptr) 0);
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


void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj)
{
        mpfr_set_d(result, omega, RND);
        mpfr_mul_d(result, result, omega_cj, RND);
}
