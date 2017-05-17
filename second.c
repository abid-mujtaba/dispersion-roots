// Calculate the second term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with s__ are internal to this module
void s__calc_coeff(mpfr_t coeff, struct Constants * const c, mpfr_t x, mpfr_t y);
void s__calc_term(mpfr_t term, struct Constants * const c, mpfr_t x, mpfr_t ic, mpfr_t * const vars);
void s__calc_term_zero(mpfr_t term, struct Constants * const c, mpfr_t * const vars);

void s__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y);

void s__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);
void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);


void calc_second(mpfr_t second, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars)
{
        s__calc_coeff(coeff, c, * vars, * (vars + 1));

        if (mpfr_cmp_ui(c->two_lambda, 0) == 0)          // Special Case - two_lambda_j == 0
                s__calc_term_zero(term, c, vars);
        else
        {
                s__calc_term(term, c, * vars, * (vars + 1), vars + 2);
        }

        mpfr_mul(second, coeff, term, RND);           // first = coeff * term
}


void s__calc_coeff(mpfr_t coeff, struct Constants * const c, mpfr_t x, mpfr_t y)
{
        mpfr_add_d(x, c->kappa, 1.5, RND);         // x = kappa + 1.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul_ui(coeff, x, 4, RND);          // coeff = 4 * gamam(c->kappa + 1.5)

        mpfr_add_d(x, c->kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul_ui(x, x, 2, RND);              // x *= 2

        mpfr_set_ui(y, 2, RND);                 // y = 2
        mpfr_sub(y, y, c->kappa, RND);             // y -= kappa
        mpfr_mul(y, y, x, RND);                 // y *= x

        mpfr_add(coeff, coeff, y, RND);         // coeff += y = 2 * gamma(c->kappa + 0.5) * (2 - kappa)

        mpfr_sub_d(x, c->kappa, 1.5, RND);         // x = kappa - 1.5
        mpfr_pow_ui(y, x, 2, RND);              // y = x^2
        mpfr_gamma(x, x, RND);                  // x = gamma(x)
        mpfr_mul(y, y, x, RND);                 // y *= x
        mpfr_mul_ui(y, y, 3, RND);              // y *= 3 = 3 * gamma(c->kappa - 1.5) * (kappa - 1.5)^2

        mpfr_set_ui(x, 1, RND);                 // x = 1
        mpfr_sub(x, x, c->kappa, RND);             // x -= kappa
        mpfr_mul(y, y, x, RND);                 // y *= x = (1 - c->kappa)

        mpfr_add(coeff, coeff, y, RND);         // coeff += y


        // Multiply with outer-most factor
        mpfr_sub_d(x, c->kappa, 1.5, RND);         // x = kappa - 1.5
        mpfr_pow_ui(x, x, 2, RND);              // x = pow(x,2)
        mpfr_mul(coeff, coeff, x, RND);         // coeff *= (c->kappa - 1.5)^2

        mpfr_mul_d(coeff, coeff, 2 * LAMBDA, RND);   // coeff *= 2 * LABDA

        mpfr_sub_d(x, c->kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_div(coeff, coeff, x, RND);         // coeff /= (c->kappa - 0.5)
}


void s__calc_term(mpfr_t term, struct Constants * const cs, mpfr_t x, mpfr_t ic, mpfr_t * const vars)
{
        struct coeffs_2f3 c;

        mpfr_set_ui(term, 1, RND);

        s__calc_coeffs_2f3_outer(& c, cs, vars);
        hyp2F3(x, c, cs->two_lambda);                     // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()


        s__calc_inner_coeff(ic, cs, * vars, * (vars + 1));


        s__calc_coeffs_2f3_inner(& c, cs, vars);
        norm_hyp2F3(x, c, cs->two_lambda);                // x = 2F3()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_sub(term, term, x, RND);                   // term -= ic * 1F2()
}


void s__calc_term_zero(mpfr_t term, struct Constants * const cs, mpfr_t * const vars)
{
        struct coeffs_2f3 c;

        if (mpfr_cmp_d(cs->kappa, 0.5) < 0)
                mpfr_printf("\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        s__calc_coeffs_2f3_outer(& c, cs, vars);

        // when two_lambda_j is zero the only surviving term is the second term of 2F3. The first term equals 1 and cancels with the 1 added to 2F3.
        // The second term has k_perp^2 and this will survive after being cancelled by 1 / k_perp^2 factor multiplied outside
        // All higher powers of k_perp^2 go to zero
        // The first term is easily calculated using the 2F3 coeffs c which create the Pochhammer symbols

        mpfr_mul(term, * c.a1, * c.a2, RND);
        mpfr_div(term, term, * c.b1, RND);
        mpfr_div(term, term, * c.b2, RND);
        mpfr_div(term, term, * c.b3, RND);

        mpfr_mul_si(term, term, -1, RND);               // 2F3 is subtracted in the term so the first term by itself should also be subtracted
}


void s__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y)
{
        if (mpfr_cmp_ui(c->omega_by_omega_c, 0) == 0)            // For om == 0 we have om / csc(cs->pi * om) = 1 / pi since limit(om -> 0) sin(om * pi) / om = pi
        {
                mpfr_set_ui(ic, 1, RND);
                mpfr_div(ic, ic, c->pi, RND);
        }
        else
        {
                mpfr_mul(ic, c->csc, c->omega_by_omega_c, RND);              // ic = csc * om
        }


        mpfr_mul_ui(ic, ic, 2, RND);    // ic *= 2
        mpfr_mul(ic, ic, c->pi, RND);      // ic *= pi

        mpfr_set_ui(x, 4, RND);                 // x = 4
        mpfr_mul_si(y, c->kappa, -1, RND);         // y = - kappa
        mpfr_pow(x, x, y, RND);                 // x = 4^(-c->kappa)
        mpfr_mul(ic, ic, x, RND);               // ic *= 4^(-c->kappa)

        mpfr_add_d(x, c->kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (c->kappa + 0.5)

        mpfr_sub_d(x, c->kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_pow(y, c->two_lambda, x, RND);      // y = two_lambda_j ^ (c->kappa - 0.5)
        mpfr_mul(ic, ic, y, RND);               // ic *= two_lambda_j ^ (c->kappa - 0.5)

        mpfr_mul_ui(x, c->kappa, 2, RND);          // x = kappa * 2
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(2 * c->kappa)

        mpfr_set_d(x, 0.5, RND);                // x = 0.5
        mpfr_sub(x, x, c->kappa, RND);             // x -= kappa
        mpfr_gamma(y, x, RND);                  // y = gamma(0.5 - c->kappa)
        mpfr_mul(ic, ic, y, RND);               // ic *= gamma(0.5 - c->kappa)

        mpfr_sub_d(x, c->kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_div(ic, ic, y, RND);               // ic /= gamma(c->kappa - 0.5)

        mpfr_add_d(x, c->kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_add(x, x, c->omega_by_omega_c, RND);                // x += om
        mpfr_gamma(y, x, RND);                  // y = gamma(x)
        mpfr_div(ic, ic, y, RND);               // ic /= gamma(c->kappa + 0.5 + c->omega_by_omega_c)
}


void s__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;
        c->b3 = vars + 4;

        mpfr_set(* c->a1, cs->kappa, RND);
        mpfr_add_d(* c->a2, cs->kappa, 1.5, RND);
        mpfr_add_d(* c->b1, cs->kappa, 0.5, RND);

        mpfr_add(* c->b2, * c->b1, cs->omega_by_omega_c, RND);
        mpfr_sub(* c->b3, * c->b1, cs->omega_by_omega_c, RND);
}


void s__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;
        c->b3 = vars + 4;

        mpfr_set_d(* c->a1, 0.5, RND);
        mpfr_set_ui(* c->a2, 2, RND);

        mpfr_set_d(* c->b1, 1.5, RND);
        mpfr_sub(* c->b1, * c->b1, cs->kappa, RND);

        mpfr_set_ui(* c->b2, 1, RND);
        mpfr_sub(* c->b2, * c->b2, cs->omega_by_omega_c, RND);

        mpfr_set_ui(* c->b3, 1, RND);
        mpfr_add(* c->b3, * c->b3, cs->omega_by_omega_c, RND) ;
}
