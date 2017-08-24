// Calculate the inner term (based on value of n) in the specie portion of the dispersion relation.

#include <mpfr.h>

#include "alpha.h"
#include "constants.h"

void term(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
        alpha(res, n, LAMBDA, c->kappa, * vars, * (vars + 1));
}
