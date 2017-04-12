// Calculate the first term

#include <mpfr.h>
#include "constants.h"


// Function prototypes. Those prefixed with f__ are internal to this module
void f__calc_coeff(mpfr_t coeff, mpfr_t kappa);
void f__calc_term(mpfr_t term, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t csc, mpfr_t pi);



void calc_first(mpfr_t first, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t csc, mpfr_t pi)
{
        mpfr_t coeff, term;
        mpfr_inits(coeff, term, (mpfr_ptr) 0);

        f__calc_coeff(coeff, kappa);
        f__calc_term(term, kappa, omega_by_omega_cj, csc, pi);

        mpfr_mul(first, coeff, term, RND);           // first = coeff * term

        mpfr_clears(coeff, term, (mpfr_ptr) 0);
}


void f__calc_coeff(mpfr_t coeff, mpfr_t kappa)
{
        // ToDo: Remove place-holder
        mpfr_set_d(coeff, 0, RND);
}


void f__calc_term(mpfr_t term, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t csc, mpfr_t pi)
{
        // ToDo: Remove place-holder
        mpfr_set_d(term, 0, RND);
}
