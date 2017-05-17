// Calculate the second term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with t__ are internal to this module
void t__calc_coeff(mpfr_t coeff, struct Constants * const c, mpfr_t x, mpfr_t y);
void t__calc_term(mpfr_t term, struct Constants * const c, mpfr_t x, mpfr_t ic, mpfr_t * const vars);
void t__calc_term_zero(mpfr_t term, struct Constants * const c, mpfr_t * const vars);

void t__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y);

void t__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);
void t__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);


void calc_third(mpfr_t third, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars)
{
        t__calc_coeff(coeff, c, * vars, * (vars + 1));

        if (mpfr_cmp_ui(c->two_lambda, 0) == 0)          // Special Case - two_lambda_j == 0
                t__calc_term_zero(term, c, vars);
        else
                t__calc_term(term, c, * vars, * (vars + 1), vars + 2);

        mpfr_mul(third, coeff, term, RND);           // third = coeff * term
}


void t__calc_coeff(mpfr_t c, struct Constants * const cs, mpfr_t x, mpfr_t y)
{
        mpfr_add(c, cs->g_k_p1_2, cs->g_k_p3_2, RND);   // c = gamma(kappa + 1/2) + gamma(kappa + 3/2)

        mpfr_mul_d(y, cs->g_k_m1_2, 0.75, RND);            // y = 0.75 * gamma(kappa _ 1/2)
        mpfr_add(c, c, y, RND);                 // c += (3/4) * gamma(kappa - 1/2)

        // Multiply with outer-most factor
        mpfr_mul_ui(c, c, 8, RND);              // c *= 8
        mpfr_mul_d(c, c, LAMBDA, RND);          // c *= LAMBDA

        mpfr_set_ui(x, 1, RND);
        mpfr_sub(x, x, cs->kappa, RND);             // x = 1 - kappa
        mpfr_mul(c, c, x, RND);                 // c * = (1 - kappa)

        mpfr_sub_d(x, cs->kappa, 1.5, RND);         // x = kappa - 3/2
        mpfr_mul(c, c, x, RND);                 // c *= (kappa - 3/2)

        mpfr_sub_d(x, cs->kappa, 0.5, RND);         // x = kappa - 1/2
        mpfr_div(c, c, x, RND);                 // c /= (kappa - 0.5)
}


void t__calc_term(mpfr_t term, struct Constants * const cs, mpfr_t x,  mpfr_t ic, mpfr_t * const vars)
{
        struct coeffs_2f3 c;

        mpfr_set_ui(term, 1, RND);

        t__calc_coeffs_2f3_outer(& c, cs, vars);
        hyp2F3(x, c, cs->two_lambda);                     // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()

        t__calc_inner_coeff(ic, cs, * vars, * (vars + 1));

        t__calc_coeffs_2f3_inner(& c, cs, vars);
        norm_hyp2F3(x, c, cs->two_lambda);                // x = 2F3()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_add(term, term, x, RND);                   // term += ic * 2F3()
}


void t__calc_term_zero(mpfr_t term, struct Constants * const cs, mpfr_t * const vars)
{
        struct coeffs_2f3 c;

        if (mpfr_cmp_d(cs->kappa, 0.5) < 0)
                mpfr_printf("\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        t__calc_coeffs_2f3_outer(& c, cs, vars);

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


void t__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y)
{
        if (mpfr_cmp_ui(c->omega_by_omega_c, 0) == 0)            // Apply L'Hoc->pital's rule to get the limit which equals 1 / pi
        {
                mpfr_set_ui(ic, 1, RND);
                mpfr_div(ic, ic, c->pi, RND);
        }
        else
        {
                mpfr_mul(ic, c->csc, c->omega_by_omega_c, RND);              // ic = csc * om
        }


        mpfr_sqrt(x, c->pi, RND);                  // x = sqrt(pi)
        mpfr_mul(ic, ic, x, RND);               // ic *= sqrt(c->pi)

        mpfr_add_d(x, c->kappa, 0.5, RND);         // x = kappa + 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (kappa + 1/2)

        mpfr_sub_d(x, c->kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_mul(ic, ic, x, RND);               // ic *= (kappa - 1/2)

        mpfr_sub_d(x, c->kappa, 1.5, RND);         // x = kappa - 3/2
        mpfr_pow(y, c->two_lambda, x, RND);      // y = two_lambda_j ^ (kappa - 3/2)
        mpfr_mul(ic, ic, y, RND);               // ic *= two_lambda_j ^ (kappa - 3/2)

        mpfr_mul(ic, ic, c->g_k_m1, RND);               // ic *= gamma(kappa - 1)

        mpfr_mul(ic, ic, c->g_5_2_mk, RND);               // ic *= gamma(5/2 - kappa)

        mpfr_sub_d(x, c->kappa, 0.5, RND);         // x = kappa - 0.5
        mpfr_add(x, x, c->omega_by_omega_c, RND);                // x += om
        mpfr_gamma(y, x, RND);                  // y = gamma(kappa - 1/2 + om)
        mpfr_mul_ui(y, y, 2, RND);              // y *= 2
        mpfr_div(ic, ic, y, RND);               // ic /= 2 * gamma(kappa - 1/2 + om)
}


void t__calc_coeffs_2f3_inner(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;
        c->b3 = vars + 4;

        mpfr_sub_ui(* c->a1, cs->kappa, 1, RND);      // a1 = kappa - 1
        mpfr_add_d(* c->a2, cs->kappa, 1.5, RND);     // a2 = kappa + 3/2

        mpfr_sub_d(* c->b1, cs->kappa, 0.5, RND);     // b1 = kappa - 1/2

        mpfr_add(* c->b2, * c->b1, cs->omega_by_omega_c, RND);        // b2 = kappa - 1/2 + om
        mpfr_sub(* c->b3, * c->b1, cs->omega_by_omega_c, RND);        // b3 = kappa - 1/2 - om
}


void t__calc_coeffs_2f3_outer(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;
        c->b3 = vars + 4;

        mpfr_set_d(* c->a1, 0.5, RND);
        mpfr_set_ui(* c->a2, 3, RND);

        mpfr_set_d(* c->b1, 2.5, RND);
        mpfr_sub(* c->b1, * c->b1, cs->kappa, RND);     // b1 = 5/2 - kappa

        mpfr_set_ui(* c->b2, 1, RND);
        mpfr_sub(* c->b2, * c->b2, cs->omega_by_omega_c, RND);        // b2 = 1 - om

        mpfr_set_ui(* c->b3, 1, RND);
        mpfr_add(* c->b3, * c->b3, cs->omega_by_omega_c, RND);        // b3 = 1 + om
}
