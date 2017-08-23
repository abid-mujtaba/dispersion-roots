#include <stdio.h>
#include <mpfr.h>

#include "alpha.h"
#include "constants.h"


void alpha1(mpfr_t result, double lambda, mpfr_t kappa);
void alpha2(mpfr_t result, double lambda, mpfr_t kappa);
void alpha3(mpfr_t result, double lambda, mpfr_t kappa);


void alpha(mpfr_t result, int n, double lambda, mpfr_t kappa)
{
    switch(n)
    {
        case 1:
            alpha1(result, lambda, kappa);
            break;

        case 2:
            alpha2(result, lambda, kappa);
            break;

        case 3:
            alpha3(result, lambda, kappa);
            break;

        default:
            fprintf(stderr, "Error: n = %d passed to alpha()", n);
    }
}


void alpha1(mpfr_t result, double lambda, mpfr_t kappa)
{
    mpfr_set_d(result, 1.234, RND);
}

void alpha2(mpfr_t result, double lambda, mpfr_t kappa)
{

}

void alpha3(mpfr_t result, double lambda, mpfr_t kappa)
{

}
