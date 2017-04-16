// Calculate the second term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with s__ are internal to this module
void s__calc_coeff(mpfr_t coeff, const mpfr_t kappa);
void s__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi);

void s__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j);

void s__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);
void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);


void calc_second(mpfr_t second, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi)
{
        mpfr_t coeff, term;
        mpfr_inits(coeff, term, (mpfr_ptr) 0);

        s__calc_coeff(coeff, kappa);
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

        mpfr_mul_d(coeff, coeff, 4 * LAMBDA, RND);   // coeff *= 4 * LABDA

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


        mpfr_set_ui(ic, 1, RND);
        // s__calc_inner_coeff(ic, csc, pi, omega_by_omega_cj, kappa, two_lambda_j);


        // s__calc_coeffs_2f3_inner(& c, kappa, omega_by_omega_cj);
        // norm_hyp1F2(x, c1f2, two_lambda_j);             // x = 1F2()
        // mpfr_mul(x, x, ic, RND);                        // x *= ic
        // mpfr_sub(term, term, x, RND);                   // term -= ic * 1F2()



        clear_coeffs_2f3(& c);
        mpfr_clears(x, ic, (mpfr_ptr) 0);
}


void s__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j)
{
        mpfr_t x, y;
        mpfr_inits(x, y, (mpfr_ptr) 0);


        mpfr_set(ic, csc, RND);                 // ic = csc

        mpfr_sqrt(x, pi, RND);                  // x = sqrt(pi)
        mpfr_mul(ic, ic, x, RND);               // ic *= sqrt(pi)
        mpfr_mul(ic, ic, om, RND);              // ic *= om

        mpfr_mul_ui(x, kappa, 2, RND);          // x = 2 * kappa
        mpfr_add_ui(x, x, 1, RND);              // x += 1
        mpfr_mul(ic, ic, x, RND);               // ic *= (2 * kappa)

        mpfr_add_d(y, kappa, 0.5, RND);         // y = kappa + 0.5
        mpfr_pow(x, two_lambda_j, y, RND);      // x = pow(two_lambda_j, y)
        mpfr_mul(ic, ic, x, RND);               // ic *= pow(two_lambda_j, kappa + 0.5)

        mpfr_add_ui(x, kappa, 1, RND);             // x = kappa + 1
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(ic, ic, x, RND);               // ic *= gamma(kappa + 1)

        mpfr_set_d(x, -0.5, RND);               // x = -0.5
        mpfr_sub(x, x, kappa, RND);                // x -= kappa
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(ic, ic, x, RND);               // ic *= gamma(-0.5 - kappa)

        mpfr_set_d(x, 1.5, RND);                // x = 1.5
        mpfr_add(x, x, kappa, RND);             // x += kappa
        mpfr_add(x, x, om, RND);                // x += omega_by_omega_cj
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul_ui(x, x, 2, RND);              // x *= 2
        mpfr_div(ic, ic, x, RND);               // ic /= 2 * gamma(kappa + 1.5 + omega_by_omega_cj)


        mpfr_clears(x, y, (mpfr_ptr) 0);
}


void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_add_ui(c->a1, kappa, 1, RND);              // a1 = kappa + 1
        mpfr_add_d(c->b1, kappa, 1.5, RND);             // b1 = kappa + 1.5

        mpfr_sub(c->b2, c->b1, om, RND);                // b2 = b1 - om  // Note: it is c.b2 which has the negative omega_by_omega_cj
        mpfr_add(c->b1, c->b1, om, RND);                // b1 += om
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
