// Functions that calculate the relevant portion of the dispersion relation in the case of kappa_j -> infinity

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function Prototypes
void calc_coeffs_2f2(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);



double specie_kappa_infinity(const double k_perp, const double omega, struct Constants * const c, mpfr_t * vars)
{
    double r;

    mpfr_t result;
    mpfr_inits(result, (mpfr_ptr) 0);

    mpfr_t * k = vars;
    mpfr_set_d(*k, 0.5, RND);       // By setting kappa = 0.5 the only term in two_lambda_j that contains kappa will become -1 (kappa - 3/2 = 1) and so its effect will be removed and a required -1 will be added

    // Calculate MPFR variables required for the three terms
    calc_omega_by_omega_cj(c->omega_by_omega_c, omega, c->omega_c);
    calc_two_lambda_j(c->two_lambda, * k, c->rho, k_perp, vars);

    // Calculate the coeffs for 2F2
    struct coeffs_2f2 c2;
    calc_coeffs_2f2(c2, c, vars + 1);

    hyp2F2(*k, c2, c->two_lambda);

    mpfr_set_d(result, 1, RND);
    mpfr_sub(result, result, *k, RND);


    // Final division
    mpfr_div(result, result, c->lambda_vc_p2, RND);

    if (k_perp != 0)            // If k_perp = 0 then the result has already been calculated by taking the limit and cancelling out the k_perp^2 term in the denominator
            mpfr_div_d(result, result, pow(k_perp, 2), RND);

    r = mpfr_get_d(result, RND);

    mpfr_clears(result, term, (mpfr_ptr) 0);
    mpfr_free_cache();              // Clear the creation of the constant pi

    return r;
}


void calc_coeffs_2f2(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;

        mpfr_set_ui(* c->a1, 1, RND);
        mpfr_set_d(* c->a2, 0.5, RND);
        mpfr_set_ui(* c->b1, 1, RND);
        mpfr_set_ui(* c->b2, 1, RND);

        mpfr_sub(* c->b1, * c->b1, cs->omega_by_omega_c, RND);
        mpfr_add(* c->b2, * c->b2, cs->omega_by_omega_c, RND);
}
