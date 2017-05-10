// Calculate the second term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with t__ are internal to this module
void t__calc_coeff(mpfr_t coeff, const mpfr_t kappa);
void t__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi);
void t__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj);

void t__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j);

void t__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);
void t__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);


void calc_third(mpfr_t third, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi)
{
        mpfr_t coeff, term;
        mpfr_inits(coeff, term, (mpfr_ptr) 0);


        t__calc_coeff(coeff, kappa);

        if (mpfr_cmp_ui(two_lambda_j, 0) == 0)          // Special Case - two_lambda_j == 0
                t__calc_term_zero(term, kappa, omega_by_omega_cj);
        else
                t__calc_term(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi);

        mpfr_mul(third, coeff, term, RND);           // third = coeff * term


        mpfr_clears(coeff, term, (mpfr_ptr) 0);
}


void t__calc_coeff(mpfr_t c, const mpfr_t kappa)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        mpfr_add_d(x, kappa, 1.5, RND);         // x = kappa + 3/2
        mpfr_gamma(c, x, RND);                  // c = gamma(kappa + 3/2)

        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 1/2
        mpfr_gamma(y, x, RND);
        mpfr_add(c, c, y, RND);                 // c += gamma(kappa + 1/2)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 1/2
        mpfr_gamma(y, x, RND);                  // y = gamma(kappa - 1/2)
        mpfr_mul_d(y, y, 0.75, RND);            // y *= 0.75
        mpfr_add(c, c, y, RND);                 // c += (3/4) * gamma(kappa - 1/2)

        // Multiply with outer-most factor
        mpfr_mul_ui(c, c, 8, RND);              // c *= 8
        mpfr_mul_d(c, c, LAMBDA, RND);          // c *= LAMBDA

        mpfr_set_ui(x, 1, RND);
        mpfr_sub(x, x, kappa, RND);             // x = 1 - kappa
        mpfr_mul(c, c, x, RND);                 // c * = (1 - kappa)

        mpfr_sub_d(x, kappa, 1.5, RND);         // x = kappa - 3/2
        mpfr_mul(c, c, x, RND);                 // c *= (kappa - 3/2)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 1/2
        mpfr_div(c, c, x, RND);                 // c /= (kappa - 0.5)

        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void t__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi)
{
        mpfr_t x, ic;
        mpfr_inits(x, ic, (mpfr_ptr) 0);
        struct coeffs_2f3 c;
        init_coeffs_2f3(& c);


        mpfr_set_ui(term, 1, RND);

        t__calc_coeffs_2f3_outer(& c, kappa, omega_by_omega_cj);
        hyp2F3(x, c, two_lambda_j);                     // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()

        t__calc_inner_coeff(ic, csc, pi, omega_by_omega_cj, kappa, two_lambda_j);

        t__calc_coeffs_2f3_inner(& c, kappa, omega_by_omega_cj);
        norm_hyp2F3(x, c, two_lambda_j);                // x = 2F3()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_add(term, term, x, RND);                   // term += ic * 2F3()


        clear_coeffs_2f3(& c);
        mpfr_clears(x, ic, (mpfr_ptr) 0);
}


void t__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t om)
{
        struct coeffs_2f3 c;
        init_coeffs_2f3(& c);

        if (mpfr_cmp_d(kappa, 0.5) < 0)
                mpfr_printf("\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        t__calc_coeffs_2f3_outer(& c, kappa, om);

        // when two_lambda_j is zero the only surviving term is the second term of 2F3. The first term equals 1 and cancels with the 1 added to 2F3.
        // The second term has k_perp^2 and this will survive after being cancelled by 1 / k_perp^2 factor multiplied outside
        // All higher powers of k_perp^2 go to zero
        // The first term is easily calculated using the 2F3 coeffs c which create the Pochhammer symbols

        mpfr_mul(term, c.a1, c.a2, RND);
        mpfr_div(term, term, c.b1, RND);
        mpfr_div(term, term, c.b2, RND);
        mpfr_div(term, term, c.b3, RND);

        mpfr_mul_si(term, term, -1, RND);               // 2F3 is subtracted in the term so the first term by itself should also be subtracted

        clear_coeffs_2f3(& c);
}


void t__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        if (mpfr_cmp_ui(om, 0) == 0)            // Apply L'Hopital's rule to get the limit which equals 1 / pi
        {
                mpfr_set_ui(ic, 1, RND);
                mpfr_div(ic, ic, pi, RND);
        }
        else
        {
                mpfr_mul(ic, csc, om, RND);              // ic = csc * om
        }


        mpfr_sqrt(x, pi, RND);                  // x = sqrt(pi)
        mpfr_mul(ic, ic, x, RND);               // ic *= sqrt(pi)

        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (kappa + 1/2)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (kappa - 1/2)

        mpfr_sub_d(x, kappa, 1.5, RND);         // x = kappa - 3/2
        mpfr_pow(y, two_lambda_j, x, RND);      // y = two_lambda_j ^ (kappa - 3/2)
        mpfr_mul(ic, ic, y, RND);               // ic *= two_lambda_j ^ (kappa - 3/2)

        mpfr_sub_ui(x, kappa, 1, RND);          // x = kappa - 1
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(kappa - 1)

        mpfr_set_d(x, 2.5, RND);                // x = 5/2
        mpfr_sub(x, x, kappa, RND);             // x -= kappa
        mpfr_gamma(y, x, RND);                  // y = gamma(5/2 - kappa)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(5/2 - kappa)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_add(x, x, om, RND);                // x += om
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_mul_ui(y, y, 2, RND);              // y *= 2
        mpfr_div(ic, ic, y, RND);               // ic /= 2 * gamma(kappa - 1/2 + om)


        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void t__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_sub_ui(c->a1, kappa, 1, RND);      // a1 = kappa - 1
        mpfr_add_d(c->a2, kappa, 1.5, RND);     // a2 = kappa + 3/2

        mpfr_sub_d(c->b1, kappa, 0.5, RND);     // b1 = kappa - 1/2

        mpfr_add(c->b2, c->b1, om, RND);        // b2 = kappa - 1/2 + om
        mpfr_sub(c->b3, c->b1, om, RND);        // b3 = kappa - 1/2 - om
}


void t__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_set_d(c->a1, 0.5, RND);
        mpfr_set_ui(c->a2, 3, RND);

        mpfr_set_d(c->b1, 2.5, RND);
        mpfr_sub(c->b1, c->b1, kappa, RND);     // b1 = 5/2 - kappa

        mpfr_set_ui(c->b2, 1, RND);
        mpfr_sub(c->b2, c->b2, om, RND);        // b2 = 1 - om

        mpfr_set_ui(c->b3, 1, RND);
        mpfr_add(c->b3, c->b3, om, RND);        // b3 = 1 + om
}
