// Calculate the inner term (based on value of n) in the specie portion of the dispersion relation.

#include <mpfr.h>
#include <stdio.h>

#include "alpha.h"
#include "constants.h"


void calculate_term(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars);


void term(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
    if (n < 1 || n > 3)
    {
        fprintf(stderr, "Error: Invalid value of n = %d (only 1, 2, 3 allowed) passed to term().", n);

        return;
    }

    calculate_term(res, n, c, * vars, * (vars + 1), vars + 2);
}


void calculate_term(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars)
{
    alpha(r, n, LAMBDA, c->kappa, * vars, * (vars + 1));
}
