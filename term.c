// Calculate the inner term (based on value of n) in the specie portion of the dispersion relation.

#include <mpfr.h>
#include <stdio.h>

#include "alpha.h"
#include "constants.h"
#include "math_utilities.h"


void second(mpfr_t r, int n, struct Constants * const c);
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars);


void term(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
    if (n < 1 || n > 3)
    {
        fprintf(stderr, "Error: Invalid value of n = %d (only 1, 2, 3 allowed) passed to term().", n);

        return;
    }

    mpfr_t * x = vars;

    // res = alpha[n] * (1 - 2F3 - third)
    mpfr_set_ui(res, 1, RND);

    second(*x, n, c);
    mpfr_sub(res, res, *x, RND);         // res = 1 - 2F3

    third(*x, n, c, * (vars + 1), * (vars + 2), vars + 3);
    mpfr_sub(res, res, *x, RND);         // res -= third

    alpha(*x, n, LAMBDA, c->kappa, * vars, * (vars + 1));
    mpfr_mul(res, res, *x, RND);         // res *= alpha[n]
}


/*
        2F3[1/2, n; n - 1/2 -k, 1 - w/w_cj, 1 + w/w_cj; 2 lambda_j]
*/
void second(mpfr_t r, int n, struct Constants * const c)
{
    // TODO: Implement
    mpfr_set_ui(r, 1, RND);
}


/*
        third = (k + 1/2)_-n / (n-1)! * sqrt(pi) csc(pi * w/w_cj)
*/
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars)
{
    mpfr_add_d(x, c->kappa, 0.5, RND);
    neg_pochammer(y, n, x, * vars);
    mpfr_mul(r, r, y, RND);                            // r *= (k + 1/2)_-n

    mpfr_div_ui(r, r, factorial(n - 1), RND);          // r /= (n - 1)!

    mpfr_mul(r, r, c->sqrt_pi, RND);            // r *= sqrt(pi)
    mpfr_mul(r, r, c->csc, RND);                // r *= csc(pi * w_j / w_ce)
}
