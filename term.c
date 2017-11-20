// Calculate the inner term (based on value of n) in the specie portion of the dispersion relation.

#include <stdio.h>
#include <mpfr.h>

#include "alpha.h"
#include "constants.h"
#include "hypergeom.h"
#include "math_utilities.h"


void term_zero(mpfr_t res, int n, struct Constants * const cs, mpfr_t * const vars);
void second(mpfr_t r, int n, struct Constants * const c, mpfr_t * vars);
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars);
void calc_second_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars);
void calc_third_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars);

// Defined in infinite_kappa.c
void term_infinite_kappa(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars);


// Calculate 1 / two_lambda_j * alpha[n] * (1 - 2F3 - third)
void term(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
    if (n < 1 || n > 3)
    {
        fprintf(stderr, "Error: Invalid value of n = %d (only 1, 2, 3 allowed) passed to term().", n);

        return;
    }

    // The special case of kappa -> infinity must be handled separately
    // TODO: Modify this to use 1 \ two_lambda_j * (1 - 2F2)
    if (mpfr_inf_p(c->kappa))
    {
        term_infinite_kappa(res, n, c, vars);

        return;
    }

    mpfr_t * x = vars;

    // if (mpfr_cmp_ui(c->two_lambda, 0) == 0)         // Special case: two_lambda_j == 0
        // term_zero(res, n, c, vars);
    // else
    // {

    // res = alpha[n] * ((1 - 2F3) - third)
    // mpfr_set_ui(res, 1, RND);

    // second(*x, n, c, vars + 1);
    // mpfr_sub(res, res, *x, RND);         // res = 1 - 2F3

    mpfr_set_ui(res, 0, RND);           // res = 0

    second(*x, n, c, vars + 1);        // 2F3 - 1
    mpfr_sub(res, res, *x, RND);        // res -= (2F3 - 1) = 1 - 2F3

    third(*x, n, c, * (vars + 1), * (vars + 2), vars + 3);
    mpfr_sub(res, res, *x, RND);         // res -= third
    // }

    alpha(*x, n, c->lambda, c->kappa, * (vars + 1), * (vars + 2));
    mpfr_mul(res, res, *x, RND);         // res *= alpha[n]
}


void term_zero(mpfr_t res, int n, struct Constants * const cs, mpfr_t * const vars)
{
    struct coeffs_2f3 c;

    if (mpfr_cmp_d(cs->kappa, 0.5) < 0)
        mpfr_fprintf(stderr, "\nError - (kappa - 1.5) < 0 - This violates the assumption used to calculate the term for k_perp = 0");

    calc_second_coeffs_2f3(& c, n, cs->kappa, cs->omega_by_omega_c, vars);

    // when two_lambda_j is zero the only surviving term is the second term of 2F3 (the third term is zero because of the two_lambda_j ^ (kappa + 3/2 - n)). The first term equals 1 and cancels with the 1 added to 2F3.
    // The second term has k_perp^2 and this will survive after being cancelled by 1 / k_perp^2 factor multiplied outside
    // All higher powers of k_perp^2 go to zero
    // The first term is easily calculated using the 2F3 coeffs c which create the Pochhammer symbols

    mpfr_mul(res, * c.a1, * c.a2, RND);
    mpfr_div(res, res, * c.b1, RND);
    mpfr_div(res, res, * c.b2, RND);
    mpfr_div(res, res, * c.b3, RND);

    // The coefficient of k_perp^2 in two_lambda_j remain after the calculation and so need to be incorporated
    mpfr_t * x = vars;
    mpfr_sub_d(*x, cs->kappa, 1.5, RND);
    mpfr_mul(res, res, *x, RND);                        // res *= (k - 3/2)
    mpfr_mul_d(res, res, 2 * cs->rho * cs->rho, RND);   // res *= 2 * rho_j^2

    mpfr_mul_si(res, res, -1, RND);               // 2F3 is subtracted in the res so the first res by itself should also be subtracted
}

/*
        1 / two_lambda_j * (1 - 2F3[1/2, n; n - 1/2 -k, 1 - w/w_cj, 1 + w/w_cj; 2 lambda_j])
*/
void second(mpfr_t r, int n, struct Constants * const c, mpfr_t * vars)
{
    struct coeffs_2f3 c_2f3;
    calc_second_coeffs_2f3(& c_2f3, n, c->kappa, c->omega_by_omega_c, vars);

    first_hyp2F3(r, c_2f3, c->two_lambda);
}

/*
        2F3[1/2, n; n - 1/2 -k, 1 - w/w_cj, 1 + w/w_cj; 2 lambda_j]
*/
void calc_second_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars)
{
    c->a1 = vars;
    c->a2 = vars + 1;
    c->b1 = vars + 2;
    c->b2 = vars + 3;
    c->b3 = vars + 4;

    mpfr_set_d(* c->a1, 0.5, RND);          // a1 = 1/2
    mpfr_set_ui(* c->a2, n, RND);           // a2 = n

    mpfr_sub_d(* c->b1, * c->a2, 0.5, RND);
    mpfr_sub(* c->b1, * c->b1, k, RND);         // b1 = n - 1/2 - k

    mpfr_set_ui(* c->b2, 1, RND);
    mpfr_sub(* c->b2, * c->b2, w, RND);         // b2 = 1 - w

    mpfr_set_ui(* c->b3, 1, RND);
    mpfr_add(* c->b3, * c->b3, w, RND);         // b3 = 1 + w
}

/*
        third = 1 / two_lambda_j * (k + 1/2)_-n / (n-1)! * sqrt(pi) csc(pi * w/w_cj) (2 lambda_j)^(k + 3/2 - n) Gamma(k + 2 - n) Gamma(n - 3/2 - k)
        / Gamma(k - n + 5/2 -w) / Gamma(k - n + 5/2 + w)
*/
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars)
{
    mpfr_add_d(x, c->kappa, 0.5, RND);
    neg_pochammer(r, n, x, * vars);             // r = (k + 1/2)_-n

    mpfr_div_ui(r, r, factorial(n - 1), RND);          // r /= (n - 1)!

    mpfr_mul(r, r, c->sqrt_pi, RND);            // r *= sqrt(pi)
    mpfr_mul(r, r, c->csc, RND);                // r *= csc(pi * w_j / w_ce)
    mpfr_mul(r, r, c->omega_by_omega_c, RND);   // r *= w_j / w_ce

    mpfr_add_ui(x, c->kappa, 2, RND);           // x = k + 2
    mpfr_sub_ui(x, x, n, RND);                  // x -= n
    Gamma(y, x);
    mpfr_mul(r, r, y, RND);                     // r *= Gamma(k + 2 - n)

    mpfr_set_ui(x, n, RND);
    mpfr_sub_d(x, x, 1.5, RND);                 // x = n - 3/2
    mpfr_sub(x, x, c->kappa, RND);              // x -= k
    Gamma(y, x);
    mpfr_mul(r, r, y, RND);                     // r *= Gamma(n - 3/2 - k)

    mpfr_sub_ui(x, c->kappa, n, RND);
    mpfr_add_d(x, x, 2.5, RND);
    mpfr_add(x, x, c->omega_by_omega_c, RND);       // x = k - n + 5/2 + w/w_ce
    Gamma(y, x);
    mpfr_div(r, r, y, RND);                     // r /= y


    mpfr_add_d(x, c->kappa, 1.5, RND);
    mpfr_sub_ui(x, x, n + 1, RND);                  // x = k + 3/2 - n - 1 (where -1 comes from the k_perp^2 in the outer denom)
    // mpfr_pow(y, c->two_lambda, x, RND);         // y = two_lambda ^ (k + 3/2 - n)
    // mpfr_mul(r, r, y, RND);                     // r *= two_lambda ^ (k + 3/2 - n)

    // Calc coeffs of inner 2F3
    struct coeffs_2f3 c_2f3;
    calc_third_coeffs_2f3(& c_2f3, n, c->kappa, c->omega_by_omega_c, vars);
    second_norm_hyp2F3(y, c_2f3, c->two_lambda, x);

    mpfr_mul(r, r, y, RND);         // r *= (two_lambda_j)^(k + 3/2 -n -1) * 2F3 / Gamma(k - n + 5/2 - w)
}


/*
        2F3[k - n + 2, k + 3/2; k + 5/2 - n; k + 5/2 - n + w; k + 5/2 -n - w; x]
*/
void calc_third_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars)
{
    c->a1 = vars;
    c->a2 = vars + 1;
    c->b1 = vars + 2;
    c->b2 = vars + 3;
    c->b3 = vars + 4;

    mpfr_sub_ui(* c->a1, k, n, RND);
    mpfr_add_ui(* c->a1, * c->a1, 2, RND);      // a1 = k - n + 2

    mpfr_add_d(* c->a2, k, 1.5, RND);           // a2 = k + 3/2

    mpfr_add_d(* c->b1, k, 2.5, RND);           // b1 = k + 5/2
    mpfr_sub_ui(* c->b1, * c->b1, n, RND);      // b1 = k + 5/2 - n

    mpfr_add(* c->b2, * c->b1, w, RND);         // b2 = k + 5/2 - n + w
    mpfr_sub(* c->b3, * c->b1, w, RND);         // b3 = k + 5/2 - n - w
}
