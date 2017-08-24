#include <stdio.h>
#include <mpfr.h>

#include "alpha.h"
#include "constants.h"


void alpha1(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y);
void alpha2(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y);
void alpha3(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y);


// ToDo: Switch to using vars instead of x and y
void alpha(mpfr_t result, int n, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{
    switch(n)
    {
        case 1:
            alpha1(result, lambda, kappa, x, y);
            break;

        case 2:
            alpha2(result, lambda, kappa, x, y);
            break;

        case 3:
            alpha3(result, lambda, kappa, x, y);
            break;

        default:
            fprintf(stderr, "Error: n = %d passed to alpha()", n);
    }
}


/*  \alpha_1 = (\kappa + 1) ( (\kappa - \frac{1}{2}) (\kappa + \frac{1}{2}) - \Lambda (\kappa - \frac{3}{2})^2 )
 */
void alpha1(mpfr_t res, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{
    mpfr_sub_d(res, kappa, 0.5, RND);       // res = (kappa - 1/2)

    mpfr_add_d(x, kappa, 0.5, RND);
    mpfr_mul(res, res, x, RND);             // res *= (kappa + 1/2)

    mpfr_sub_d(y, kappa, 1.5, RND);         // y = (kappa - 3/2)
    mpfr_mul(y, y, x, RND);                 // y *= (kappa - 3/2)  -  squared
    mpfr_mul_d(y, y, lambda, RND);          // y *= Lambda

    mpfr_sub(res, res, y, RND);             // res -= y

    mpfr_add_ui(x, kappa, 1, RND);
    mpfr_mul(res, res, x, RND);             // res *= (kappa + 1)
}

void alpha2(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{

}

void alpha3(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{

}
