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
void specie(mpfr_t result, double k_perp, double omega, struct Constants * const cj, mpfr_t * vars);

// The following functions are declared here but are defined elsewhere
void term(mpfr_t result, int n, struct Constants * const c, mpfr_t * const vars);
void calc_first(mpfr_t first, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_second(mpfr_t second, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);
void calc_third(mpfr_t third, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars);


double D(const double k_perp, const double omega)
{
        int p = 0;
        double r;

        mpfr_t result, x;
        mpfr_inits(result, x, (mpfr_ptr) 0);

        // We start by setting the default precision for MPFR variables based on the value of k_perp. The larger it is the higher the precision required.
        p = 1 + (int) (k_perp / DOUBLE_PRECISION_DELTA);

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
        mpfr_set_ui(result, 1, RND);
        specie(x, k_perp, omega, &cc, vars);
        mpfr_add(result, result, x, RND);

        specie(x, k_perp, omega, &ch, vars);
        mpfr_add(result, result, x, RND);

        r = mpfr_get_d(result, RND);

        // Clear the variables
        for (int i = 0; i < NUM_MPFR_VARIABLES; ++i)
                mpfr_clear(vars[i]);

        // Clear variables inside the constants
        clear_constants(& cc);
        clear_constants(& ch);

        mpfr_clears(result, x, (mpfr_ptr) 0);
        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)

        if (isnan(r))
            fprintf(stderr, "\nWarning - D(%.2f, %.2f) = NAN", k_perp, omega);

        return r;
}


void specie(mpfr_t result, const double k_perp, const double omega, struct Constants * const c, mpfr_t * vars)
{
        // TODO: Handle the case of infinite kappa

        mpfr_set_ui(result, 0, RND);                           // res = 0

        // mpfr_t term;
        // mpfr_inits(term, (mpfr_ptr) 0);

        // Pre-Calculate MPFR variables required for internal calculation
        calc_omega_by_omega_cj(c->omega_by_omega_c, omega, c->omega_c);
        calc_two_lambda_j(c->two_lambda, c->kappa, c->rho, k_perp, vars);

        mpfr_mul(c->csc, c->omega_by_omega_c, c->pi, RND);    // csc = omega_by_omega_cj * pi
        mpfr_csc(c->csc, c->csc, RND);                        // csc = cosec( csc )

        mpfr_t *t = vars;          // Temporary variable

        for (int n = 1; n <=3; ++n)
        {
            term(*t, n, c, vars + 1);       // Calculate term for specified value of n
            mpfr_add(result, result, *t, RND);        // result += term[n]
        }


        // Final division
        mpfr_div(result, result, c->lambda_vc_p2, RND);

        // TODO: Ensure that this is functioning correctly in the consolidated approach (using term())
        if (k_perp != 0)            // If k_perp = 0 then the result has already been calculated by taking the limit and cancelling out the k_perp^2 term in the denominator
                mpfr_div_d(result, result, pow(k_perp, 2), RND);
}


// two_lambda_j = 2 (k - 3/2) (k_perp * rho_j)^2
void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp, mpfr_t * vars)
{
        mpfr_t * x = vars;              // x equals the first pointer in vars, which has already been initialized

        mpfr_set_d(result, k_perp, RND);                // result = k_perp;
        mpfr_mul_d(result, result, rho_j, RND);         // result *= rho_j;
        mpfr_pow_ui(result, result, 2, RND);            // result = pow(result, 2);
        mpfr_mul_ui(result, result, 2, RND);            // result *= 2

        // We de-reference the pointer x to gain access to the mpfr_t variable it points to
        // If kappa -> infinity it is removed fom two_lambda_j since it will cancel out with the kappa term in the 2F3 making it a 2F2 but the cancellation leaves a factor of -1 which we incorporate here
        if (mpfr_inf_p(kappa_j))
            mpfr_mul_si(result, result, -1, RND);       // result *= -1
        else
        {
            mpfr_sub_d(*x, kappa_j, 1.5, RND);
            mpfr_mul(result, result, *x, RND);           // result *= (kappa - 1.5)
        }
}


void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj)
{
        mpfr_set_d(result, omega, RND);
        mpfr_mul_d(result, result, omega_cj, RND);
}
