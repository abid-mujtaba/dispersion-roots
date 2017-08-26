// Functions that calculate the relevant portion of the dispersion relation in the case of kappa_j -> infinity

#include <stdio.h>
#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


void calc_coeffs_2f2(struct coeffs_2f2 * const c, int n, const mpfr_t w, mpfr_t * vars);
void alpha_infinite_kappa(mpfr_t res, int n, double lambda);


// In the case of infinite kappa only the first 2F3 survives as a 2F2 (without the kappa term)
void term_infinite_kappa(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
    mpfr_t * x = vars;

    // Calculate coeffs and value of 2F2
    struct coeffs_2f2 cf;

    calc_coeffs_2f2(& cf, n, c->omega_by_omega_c, vars);

    // Deal with the special case of k_perp = 0
    if (mpfr_cmp_ui(c->two_lambda, 0) == 0)         // Special case: two_lambda_j == 0
    {
        // In this case only the first term of the 2F2 survives with a negative sign
        mpfr_mul(res, * cf.a1, * cf.a2, RND);
        mpfr_div(res, res, * cf.b1, RND);
        mpfr_div(res, res, * cf.b2, RND);
        mpfr_mul_si(res, res, -1, RND);
    }
    else
    {
        mpfr_set_ui(res, 1, RND);
        hyp2F2(*x, cf, c->two_lambda);
        mpfr_sub(res, res, *x, RND);
    }

    // Multiply with alpha
    alpha_infinite_kappa(*x, n, LAMBDA);
    mpfr_mul(res, res, *x, RND);
}


void alpha_infinite_kappa(mpfr_t res, int n, double lambda)
{
    switch (n)
    {
        case 1:
            mpfr_set_ui(res, 1, RND);
            mpfr_sub_d(res, res, lambda, RND);
            break;

        case 2:
            mpfr_set_d(res, 2 * lambda, RND);
            break;

        case 3:
            mpfr_set_d(res, 8 * lambda, RND);
            break;

        default:
            fprintf(stderr, "\nError - Invalid value of n = %d while calculating term for kappa -> inf.", n);
    }
}


/*
        2F2[1/2, n; 1 - w/w_cj, 1 + w/w_cj; 2 lambda_j]
*/
void calc_coeffs_2f2(struct coeffs_2f2 * const c, int n, const mpfr_t w, mpfr_t * vars)
{
    c->a1 = vars;
    c->a2 = vars + 1;
    c->b1 = vars + 2;
    c->b2 = vars + 3;

    mpfr_set_d(* c->a1, 0.5, RND);          // a1 = 1/2
    mpfr_set_ui(* c->a2, n, RND);           // a2 = n

    mpfr_set_ui(* c->b1, 1, RND);
    mpfr_sub(* c->b1, * c->b1, w, RND);         // b1 = 1 - w

    mpfr_set_ui(* c->b2, 1, RND);
    mpfr_add(* c->b2, * c->b2, w, RND);         // b2 = 1 + w
}
