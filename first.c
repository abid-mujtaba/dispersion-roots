// Calculate the first term

#include <stdio.h>
#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with f__ are internal to this module
void f__calc_coeff(mpfr_t coeff, struct Constants * const c, mpfr_t x, mpfr_t y);
void f__calc_term(mpfr_t term, struct Constants * const c, mpfr_t x, mpfr_t ic, mpfr_t * const vars);
void f__calc_term_zero(mpfr_t term, struct Constants * const c, mpfr_t * const vars);

void f__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y);

void f__calc_coeffs_1f2(struct coeffs_1f2 * const c, struct Constants * const cs, mpfr_t * const vars);
void f__calc_coeffs_2f3(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars);


void calc_first(mpfr_t first, struct Constants * const c, mpfr_t coeff, mpfr_t term, mpfr_t * const vars)
{
        f__calc_coeff(coeff, c, * vars, * (vars + 1));

        if (mpfr_cmp_ui(c->two_lambda, 0) == 0)          // Special Case - two_lambda_j == 0
                f__calc_term_zero(term, c, vars);
        else
                f__calc_term(term, c, * vars, * (vars + 1), vars + 2);

        mpfr_mul(first, coeff, term, RND);           // first = coeff * term
}


void f__calc_coeff(mpfr_t coeff, struct Constants * const c, mpfr_t x, mpfr_t y)
{
        mpfr_sub_d(y, c->kappa, 1.5, RND);         // y = kappa - 1.5
        mpfr_mul(y, y, y, RND);                 // y = (kappa - 1.5)^2
        mpfr_mul_d(y, y, LAMBDA, RND);          // y *= LAMBDA

        mpfr_add_ui(coeff, c->kappa, 1, RND);      // coeff = kappa + 1
        mpfr_mul(coeff, coeff, c->g_k_p3_2, RND);         // coeff *= x

        mpfr_mul(x, c->g_k_p1_2, y, RND);                 // x = y * gamma(kappa + 1/2)
        mpfr_mul_ui(x, x, 4, RND);              // x *= 4

        mpfr_sub(coeff, coeff, x, RND);         // coeff -= x

        mpfr_mul(x, c->g_k_m1_2, y, RND);       // x = y * gamma(kappa - 1/2)
        mpfr_mul_ui(x, x, 3, RND);              // x *= 3

        mpfr_set_d(y, 1, RND);                  // y = 1        (new def)
        mpfr_sub(y, y, c->kappa, RND);             // y -= kappa

        mpfr_mul(x, x, y, RND);                 // x *= y

        mpfr_sub(coeff, coeff, x, RND);         // coeff -= x
}


void f__calc_term(mpfr_t term, struct Constants * const c, mpfr_t x, mpfr_t ic, mpfr_t * const vars)
{
        mpfr_set_ui(term, 1, RND);

        f__calc_inner_coeff(ic, c, * vars, * (vars + 1));

        struct coeffs_1f2 c1f2;
        struct coeffs_2f3 c2f3;

        f__calc_coeffs_1f2(& c1f2, c, vars);
        norm_hyp1F2(x, c1f2, c->two_lambda);             // x = 1F2()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_sub(term, term, x, RND);                   // term -= ic * 1F2()

        f__calc_coeffs_2f3(& c2f3, c, vars);
        hyp2F3(x, c2f3, c->two_lambda);                  // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()
}


void f__calc_term_zero(mpfr_t term, struct Constants * const cs, mpfr_t * const vars)
{
        struct coeffs_2f3 c;

        if (mpfr_cmp_d(cs->kappa, 0.5) < 0)
                mpfr_fprintf(stderr, "\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        f__calc_coeffs_2f3(& c, cs, vars);

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


void f__calc_inner_coeff(mpfr_t ic, struct Constants * const c, mpfr_t x, mpfr_t y)
{
        if (mpfr_cmp_ui(c->omega_by_omega_c, 0) == 0)            // Apply L'Hopital's rule to get the limit which equals 1 / pi
        {
                mpfr_set_ui(ic, 1, RND);
                mpfr_div(ic, ic, c->pi, RND);
        }
        else
        {
                mpfr_mul(ic, c->csc, c->omega_by_omega_c, RND);              // ic = csc * om
        }

        mpfr_sqrt(x, c->pi, RND);                  // x = sqrt(pi)
        mpfr_mul(ic, ic, x, RND);               // ic *= sqrt(pi)

        mpfr_mul_ui(x, c->kappa, 2, RND);          // x = 2 * kappa
        mpfr_add_ui(x, x, 1, RND);              // x += 1
        mpfr_mul(ic, ic, x, RND);               // ic *= (2 * kappa)

        mpfr_add_d(y, c->kappa, 0.5, RND);         // y = kappa + 0.5
        mpfr_pow(x, c->two_lambda, y, RND);      // x = pow(two_lambda_j, y)
        mpfr_mul(ic, ic, x, RND);               // ic *= pow(two_lambda_j, kappa + 0.5)

        mpfr_mul(ic, ic, c->g_k_p1, RND);           // ic *= gamma(kappa + 1)
        mpfr_mul(ic, ic, c->g_m1_2_mk, RND);        // ic *= gamma(-1/2 - kappa)

        mpfr_set_d(x, 1.5, RND);                // x = 1.5
        mpfr_add(x, x, c->kappa, RND);             // x += kappa
        mpfr_add(x, x, c->omega_by_omega_c, RND);                // x += omega_by_omega_cs
        mpfr_gamma(x, x, RND);                  // x = gamma(kappa + om + 3/2)
        mpfr_mul_ui(x, x, 2, RND);              // x *= 2
        mpfr_div(ic, ic, x, RND);               // ic /= 2 * gamma(kappa + 1.5 + omega_by_omega_cs)
}


void f__calc_coeffs_1f2(struct coeffs_1f2 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->b1 = vars + 1;
        c->b2 = vars + 2;

        mpfr_add_ui(* c->a1, cs->kappa, 1, RND);              // a1 = kappa + 1
        mpfr_add_d(* c->b1, cs->kappa, 1.5, RND);             // b1 = kappa + 1.5

        mpfr_sub(* c->b2, * c->b1, cs->omega_by_omega_c, RND);                // b2 = b1 - om  // Note: it is c.b2 which has the negative omega_by_omega_cs
        mpfr_add(* c->b1, * c->b1, cs->omega_by_omega_c, RND);                // b1 += om
}


void f__calc_coeffs_2f3(struct coeffs_2f3 * const c, struct Constants * const cs, mpfr_t * const vars)
{
        c->a1 = vars;
        c->a2 = vars + 1;
        c->b1 = vars + 2;
        c->b2 = vars + 3;
        c->b3 = vars + 4;

        mpfr_set_ui(* c->a1, 1, RND);
        mpfr_set_d(* c->a2, 0.5, RND);
        mpfr_set_d(* c->b1, 0.5, RND);
        mpfr_set_ui(* c->b2, 1, RND);

        mpfr_sub(* c->b1, * c->b1, cs->kappa, RND);
        mpfr_sub(* c->b3, * c->b2, cs->omega_by_omega_c, RND);
        mpfr_add(* c->b2, * c->b2, cs->omega_by_omega_c, RND);
}
