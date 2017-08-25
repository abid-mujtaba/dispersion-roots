// Calculate the inner term (based on value of n) in the specie portion of the dispersion relation.

#include <mpfr.h>
#include <stdio.h>

#include "alpha.h"
#include "constants.h"
#include "hypergeom.h"
#include "math_utilities.h"


void second(mpfr_t r, int n, struct Constants * const c, mpfr_t * vars);
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars);
void calc_second_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars);
void calc_third_coeffs_2f3(struct coeffs_2f3 * const c, int n, const mpfr_t k, const mpfr_t w, mpfr_t * vars);


void term(mpfr_t res, int n, struct Constants * const c, mpfr_t * const vars)
{
    if (n < 1 || n > 3)
    {
        fprintf(stderr, "Error: Invalid value of n = %d (only 1, 2, 3 allowed) passed to term().", n);

        return;
    }

    mpfr_t * x = vars;

    // res = alpha[n] * (1 - 2F3 - third)
    mpfr_set_ui(res, 1, RND);

    second(*x, n, c, vars + 1);
    mpfr_sub(res, res, *x, RND);         // res = 1 - 2F3

    third(*x, n, c, * (vars + 1), * (vars + 2), vars + 3);
    mpfr_sub(res, res, *x, RND);         // res -= third

    alpha(*x, n, LAMBDA, c->kappa, * (vars + 1), * (vars + 2));
    mpfr_mul(res, res, *x, RND);         // res *= alpha[n]
}


/*
        2F3[1/2, n; n - 1/2 -k, 1 - w/w_cj, 1 + w/w_cj; 2 lambda_j]
*/
void second(mpfr_t r, int n, struct Constants * const c, mpfr_t * vars)
{
    struct coeffs_2f3 c_2f3;
    calc_second_coeffs_2f3(& c_2f3, n, c->kappa, c->omega_by_omega_c, vars);

    norm_hyp2F3(r, c_2f3, c->two_lambda);
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
        third = (k + 1/2)_-n / (n-1)! * sqrt(pi) csc(pi * w/w_cj) (2 lambda_j)^(k + 3/2 - n) Gamma(k + 2 - n) Gamma(n - 3/2 - k)
        / Gamma(k - n + 5/2 -w) / Gamma(k - n + 5/2 + w)
*/
void third(mpfr_t r, int n, struct Constants * const c, mpfr_t x, mpfr_t y, mpfr_t * vars)
{
    mpfr_add_d(x, c->kappa, 0.5, RND);
    neg_pochammer(y, n, x, * vars);
    mpfr_mul(r, r, y, RND);                            // r *= (k + 1/2)_-n

    mpfr_div_ui(r, r, factorial(n - 1), RND);          // r /= (n - 1)!

    mpfr_mul(r, r, c->sqrt_pi, RND);            // r *= sqrt(pi)
    mpfr_mul(r, r, c->csc, RND);                // r *= csc(pi * w_j / w_ce)

    mpfr_add_d(x, c->kappa, 1.5, RND);
    mpfr_sub_ui(x, x, n, RND);                  // x = k + 3/2 - n
    mpfr_pow(y, c->two_lambda, x, RND);         // y = two_lambda ^ (k + 3/2 - n)
    mpfr_mul(r, r, y, RND);                     // r *= two_lambda ^ (k + 3/2 - n)

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

    // Calc coeffs of inner 2F3
    struct coeffs_2f3 c_2f3;
    calc_third_coeffs_2f3(& c_2f3, n, c->kappa, c->omega_by_omega_c, vars);
    norm_hyp2F3(x, c_2f3, c->two_lambda);

    mpfr_mul(r, r, x, RND);         // r *= 2F3 / Gamma(k - n + 5/2 - w)
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
