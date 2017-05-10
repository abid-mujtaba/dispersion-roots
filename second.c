// Calculate the second term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with s__ are internal to this module
void s__calc_coeff(mpfr_t coeff, const mpfr_t kappa);
void s__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi);
void s__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj);

void s__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j);

void s__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);
void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);


void calc_second(mpfr_t second, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi)
{
        mpfr_t coeff, term;
        mpfr_inits(coeff, term, (mpfr_ptr) 0);

        s__calc_coeff(coeff, kappa);

        if (mpfr_cmp_ui(two_lambda_j, 0) == 0)          // Special Case - two_lambda_j == 0
                s__calc_term_zero(term, kappa, omega_by_omega_cj);
        else
                s__calc_term(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi);

        mpfr_mul(second, coeff, term, RND);           // first = coeff * term

        mpfr_clears(coeff, term, (mpfr_ptr) 0);
}


void s__calc_coeff(mpfr_t coeff, const mpfr_t kappa)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        mpfr_add_d(x, kappa, 1.5, RND);         // x = kappa + 1.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul_ui(coeff, x, 4, RND);          // coeff = 4 * gamam(kappa + 1.5)

        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul_ui(x, x, 2, RND);              // x *= 2

        mpfr_set_ui(y, 2, RND);                 // y = 2
        mpfr_sub(y, y, kappa, RND);             // y -= kappa
        mpfr_mul(y, y, x, RND);                 // y *= x

        mpfr_add(coeff, coeff, y, RND);         // coeff += y = 2 * gamma(kappa + 0.5) * (2 - kappa)

        mpfr_sub_d(x, kappa, 1.5, RND);         // x = kappa - 1.5
        mpfr_pow_ui(y, x, 2, RND);              // y = x^2
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(y, y, x, RND);                 // y *= x
        mpfr_mul_ui(y, y, 3, RND);              // y *= 3 = 3 * gamma(kappa - 1.5) * (kappa - 1.5)^2

        mpfr_set_ui(x, 1, RND);                 // x = 1
        mpfr_sub(x, x, kappa, RND);             // x -= kappa
        mpfr_mul(y, y, x, RND);                 // y *= x = (1 - kappa)

        mpfr_add(coeff, coeff, y, RND);         // coeff += y


        // Multiply with outer-most factor
        mpfr_sub_d(x, kappa, 1.5, RND);         // x = kappa - 1.5
        mpfr_pow_ui(x, x, 2, RND);              // x = pow(x,2)
        mpfr_mul(coeff, coeff, x, RND);         // coeff *= (kappa - 1.5)^2

        mpfr_mul_d(coeff, coeff, 2 * LAMBDA, RND);   // coeff *= 2 * LABDA

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_div(coeff, coeff, x, RND);         // coeff /= (kappa - 0.5)


        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void s__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi)
{
        mpfr_t x, ic;
        mpfr_inits(x, ic, (mpfr_ptr) 0);

        struct coeffs_2f3 c;
        init_coeffs_2f3(& c);

        mpfr_set_ui(term, 1, RND);

        s__calc_coeffs_2f3_outer(& c, kappa, omega_by_omega_cj);
        hyp2F3(x, c, two_lambda_j);                     // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()


        s__calc_inner_coeff(ic, csc, pi, omega_by_omega_cj, kappa, two_lambda_j);


        s__calc_coeffs_2f3_inner(& c, kappa, omega_by_omega_cj);
        norm_hyp2F3(x, c, two_lambda_j);                // x = 2F3()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_sub(term, term, x, RND);                   // term -= ic * 1F2()



        clear_coeffs_2f3(& c);
        mpfr_clears(x, ic, (mpfr_ptr) 0);
}


void s__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t om)
{
        struct coeffs_2f3 c;
        init_coeffs_2f3(& c);

        if (mpfr_cmp_d(kappa, 0.5) < 0)
                mpfr_printf("\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        s__calc_coeffs_2f3_outer(& c, kappa, om);

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


void s__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        mpfr_mul(ic, csc, om, RND);     // ic = csc * om
        mpfr_mul_ui(ic, ic, 2, RND);    // ic *= 2
        mpfr_mul(ic, ic, pi, RND);      // ic *= pi

        mpfr_set_ui(x, 4, RND);                 // x = 4
        mpfr_mul_si(y, kappa, -1, RND);         // y = - kappa
        mpfr_pow(x, x, y, RND);                 // x = 4^(-kappa)
        mpfr_mul(ic, ic, x, RND);               // ic *= 4^(-kappa)

        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (kappa + 0.5)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_pow(y, two_lambda_j, x, RND);      // y = two_lambda_j ^ (kappa - 0.5)
        mpfr_mul(ic, ic, y, RND);               // ic *= two_lambda_j ^ (kappa - 0.5)

        mpfr_mul_ui(x, kappa, 2, RND);          // x = kappa * 2
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(2 * kappa)

        mpfr_set_d(x, 0.5, RND);                // x = 0.5
        mpfr_sub(x, x, kappa, RND);             // x -= kappa
        mpfr_gamma(y, x, RND);                  // y = gamma(0.5 - kappa)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(0.5 - kappa)

        mpfr_sub_d(x, kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_div(ic, ic, y, RND);               // ic /= gamma(kappa - 0.5)

        mpfr_add_d(x, kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_add(x, x, om, RND);                // x += om
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_div(ic, ic, y, RND);               // ic /= gamma(kappa + 0.5 + om)


        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_set(c->a1, kappa, RND);
        mpfr_add_d(c->a2, kappa, 1.5, RND);
        mpfr_add_d(c->b1, kappa, 0.5, RND);

        mpfr_add(c->b2, c->b1, om, RND);
        mpfr_sub(c->b3, c->b1, om, RND);
}


void s__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_set_d(c->a1, 0.5, RND);
        mpfr_set_ui(c->a2, 2, RND);

        mpfr_set_d(c->b1, 1.5, RND);
        mpfr_sub(c->b1, c->b1, kappa, RND);

        mpfr_set_ui(c->b2, 1, RND);
        mpfr_sub(c->b2, c->b2, om, RND);

        mpfr_set_ui(c->b3, 1, RND);
        mpfr_add(c->b3, c->b3, om, RND) ;
}
