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


/*  alpha_1 = (k + 1) ( (k - 1/2) (k + 1/2) - \Lambda (k - 3/2)^2 )
 */
void alpha1(mpfr_t res, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{
    mpfr_sub_d(res, kappa, 0.5, RND);       // res = (k - 1/2)

    mpfr_add_d(x, kappa, 0.5, RND);
    mpfr_mul(res, res, x, RND);             // res *= (k + 1/2)

    mpfr_sub_d(y, kappa, 1.5, RND);         // y = (k - 3/2)
    mpfr_mul(y, y, x, RND);                 // y *= (k - 3/2)  -  squared
    mpfr_mul_d(y, y, lambda, RND);          // y *= Lambda

    mpfr_sub(res, res, y, RND);             // res -= y

    mpfr_add_ui(x, kappa, 1, RND);
    mpfr_mul(res, res, x, RND);             // res *= (k + 1)
}


/* alpha_2 = -2 Lambda (k - 3/2)^2 / (k - 1/2) * ( 2 (k - 1/2) (k + 1) - 3 (k - 1) (k - 3/2))
*/
void alpha2(mpfr_t res, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{
    mpfr_sub_d(res, kappa, 0.5, RND);       // res = (k - 1/2)
    mpfr_mul_ui(res, res, 2, RND);          // res *= 2

    mpfr_add_ui(x, kappa, 3, RND);
    mpfr_mul(res, res, x, RND);             // res *= (k + 3)

    mpfr_sub_ui(y, kappa, 1, RND);          // y = (k - 1)
    mpfr_mul_ui(y, y, 3, RND);              // y *= 3

    mpfr_sub_d(x, kappa, 1.5, RND);
    mpfr_mul(y, y, x, RND);                 // y *= (k - 3/2)

    mpfr_sub(res, res, y, RND);             // res -= y
}

void alpha3(mpfr_t result, double lambda, mpfr_t kappa, mpfr_t x, mpfr_t y)
{

}
