// Calculate the first term

#include <mpfr.h>
#include "constants.h"


// Function prototypes. Those prefixed with f__ are internal to this module
void f__calc_coeff(mpfr_t coeff, mpfr_t kappa);
void f__calc_term(mpfr_t term, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi);
void f__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j);



void calc_first(mpfr_t first, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi)
{
        mpfr_t coeff, term;
        mpfr_inits(coeff, term, (mpfr_ptr) 0);

        f__calc_coeff(coeff, kappa);
        f__calc_term(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi);

        mpfr_mul(first, coeff, term, RND);           // first = coeff * term

        mpfr_clears(coeff, term, (mpfr_ptr) 0);
}


void f__calc_coeff(mpfr_t coeff, mpfr_t kappa)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        mpfr_sub_d(y, kappa, 1.5, RND);         // y = kappa - 1.5
        mpfr_mul(y, y, y, RND);                 // y = (kappa - 1.5)^2
        mpfr_mul_d(y, y, LAMBDA, RND);          // y *= LAMBDA


        mpfr_add_ui(coeff, kappa, 1, RND);      // coeff = kappa + 1

        mpfr_add_d(x, kappa, 1.5, RND);         // x = kappa + 1.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)

        mpfr_mul(coeff, coeff, x, RND);         // coeff *= x


        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(x, x, y, RND);                 // x *= y
        mpfr_mul_ui(x, x, 4, RND);              // x *= 4

        mpfr_sub(coeff, coeff, x, RND);         // coeff -= x


        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(x, x, y, RND);                 // x *= y
        mpfr_mul_ui(x, x, 3, RND);              // x *= 3

        mpfr_set_d(y, 1, RND);                  // y = 1        (new def)
        mpfr_sub(y, y, kappa, RND);             // y -= kappa

        mpfr_mul(x, x, y, RND);                 // x *= y

        mpfr_sub(coeff, coeff, x, RND);         // coeff -= x


        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void f__calc_term(mpfr_t term, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi)
{
        mpfr_t x, ic;
        mpfr_inits(x, ic, (mpfr_ptr) 0);

        f__calc_inner_coeff(ic, csc, pi, omega_by_omega_cj, kappa, two_lambda_j);

        // ToDo: Remove place-holder
        mpfr_mul_d(term, ic, 0, RND);


        mpfr_clears(x, ic, (mpfr_ptr) 0);
}


void f__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j)
{
        // ToDo: Remove place-holder
        mpfr_set_d(ic, 0, RND);
}
