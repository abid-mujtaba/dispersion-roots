// Calculate the first term

#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Function prototypes. Those prefixed with f__ are internal to this module
void f__calc_coeff(mpfr_t coeff, const mpfr_t kappa, mpfr_t x, mpfr_t y);
void f__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi, mpfr_t x, mpfr_t ic, mpfr_t * const vars);
void f__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj);

void f__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j, mpfr_t x, mpfr_t y);

void f__calc_coeffs_1f2(struct coeffs_1f2 * const c, const mpfr_t kappa, const mpfr_t om, mpfr_t * const vars);
void f__calc_coeffs_2f3(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om);


void calc_first(mpfr_t first, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi, mpfr_t coeff, mpfr_t term, mpfr_t * const vars)
{
        f__calc_coeff(coeff, kappa, * vars, * (vars + 1));

        if (mpfr_cmp_ui(two_lambda_j, 0) == 0)          // Special Case - two_lambda_j == 0
                f__calc_term_zero(term, kappa, omega_by_omega_cj);
        else
                f__calc_term(term, kappa, omega_by_omega_cj, two_lambda_j, csc, pi, * vars, * (vars + 1), vars + 2);

        mpfr_mul(first, coeff, term, RND);           // first = coeff * term
}


void f__calc_coeff(mpfr_t coeff, const mpfr_t kappa, mpfr_t x, mpfr_t y)
{
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
}


void f__calc_term(mpfr_t term, const mpfr_t kappa, const mpfr_t omega_by_omega_cj, const mpfr_t two_lambda_j, const mpfr_t csc, const mpfr_t pi, mpfr_t x, mpfr_t ic, mpfr_t * const vars)
{
        mpfr_set_ui(term, 1, RND);

        f__calc_inner_coeff(ic, csc, pi, omega_by_omega_cj, kappa, two_lambda_j, * vars, * (vars + 1));

        struct coeffs_1f2 c1f2;
        struct coeffs_2f3 c2f3;

        f__calc_coeffs_1f2(& c1f2, kappa, omega_by_omega_cj, vars);
        norm_hyp1F2(x, c1f2, two_lambda_j);             // x = 1F2()
        mpfr_mul(x, x, ic, RND);                        // x *= ic
        mpfr_sub(term, term, x, RND);                   // term -= ic * 1F2()

        init_coeffs_2f3(& c2f3);
        f__calc_coeffs_2f3(& c2f3, kappa, omega_by_omega_cj);
        hyp2F3(x, c2f3, two_lambda_j);                  // x = 2F3()
        mpfr_sub(term, term, x, RND);                   // term -= 2F3()


        clear_coeffs_2f3(& c2f3);
}


void f__calc_term_zero(mpfr_t term, const mpfr_t kappa, const mpfr_t om)
{
        struct coeffs_2f3 c;
        init_coeffs_2f3(& c);

        if (mpfr_cmp_d(kappa, 0.5) < 0)
                mpfr_printf("\nWarning - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

        f__calc_coeffs_2f3(& c, kappa, om);

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


void f__calc_inner_coeff(mpfr_t ic, const mpfr_t csc, const mpfr_t pi, const mpfr_t om, const mpfr_t kappa, const mpfr_t two_lambda_j, mpfr_t x, mpfr_t y)
{
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
}


void f__calc_coeffs_1f2(struct coeffs_1f2 * const c, const mpfr_t kappa, const mpfr_t om, mpfr_t * const vars)
{
        c->a1 = vars;
        c->b1 = vars + 1;
        c->b2 = vars + 2;

        mpfr_add_ui(* c->a1, kappa, 1, RND);              // a1 = kappa + 1
        mpfr_add_d(* c->b1, kappa, 1.5, RND);             // b1 = kappa + 1.5

        mpfr_sub(* c->b2, * c->b1, om, RND);                // b2 = b1 - om  // Note: it is c.b2 which has the negative omega_by_omega_cj
        mpfr_add(* c->b1, * c->b1, om, RND);                // b1 += om
}


void f__calc_coeffs_2f3(struct coeffs_2f3 * const c, const mpfr_t kappa, const mpfr_t om)
{
        mpfr_set_ui(c->a1, 1, RND);
        mpfr_set_d(c->a2, 0.5, RND);
        mpfr_set_d(c->b1, 0.5, RND);
        mpfr_set_ui(c->b2, 1, RND);

        mpfr_sub(c->b1, c->b1, kappa, RND);
        mpfr_sub(c->b3, c->b2, om, RND);
        mpfr_add(c->b2, c->b2, om, RND);
}
