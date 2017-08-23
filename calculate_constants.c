/*
 * Calculate the constants derived from the ones defined in constants.h and store
 * them in derived.h.
 *
 * This process will be carried out before compilation so that these values are
 * pre-compiled in to the executable and the calculations won't have to be
 * carried out over and over again.
 */

 // ToDo: Remove caching of unneeded gamma functions
 //       Use recursion where necessary to reduce calculations

#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "roots.h"
#include <mpfr.h>


// Store the mpfr var x with name in the output file:  #define <name> "<value of x>"
void foutput(FILE * fout, mpfr_t x, char * name);

void calc_lambda_vcj_p2(mpfr_t result, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x);
void calc_gamma(FILE * fout, mpfr_t res, mpfr_t kappa, double delta, char * name);
void calc_gamma_minus(FILE * fout, mpfr_t res, mpfr_t kappa, double delta, char * name);


int main(void)
{
    // Calculate the max needed precision from max value of k_perp
    // All constants are calculated with this max precision since downgrading is easy
    int p = 1 + (int) (K_PERP_MAX / 30);
    mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

    mpfr_t res, kappa, x;
    mpfr_inits(res, kappa, x, (mpfr_t *) 0);

    FILE * fout = fopen("derived.h", "w");


    // Hot specie calculations
    mpfr_set_d(kappa, KAPPA_H, RND);

    calc_gamma(fout, res, kappa, 0.5, "G_K_P1_2_H");
    calc_gamma(fout, res, kappa, -0.5, "G_K_M1_2_H");
    calc_gamma(fout, res, kappa, 1, "G_K_P1_H");
    calc_gamma(fout, res, kappa, -1, "G_K_M1_H");
    calc_gamma(fout, res, kappa, 1.5, "G_K_P3_2_H");
    calc_gamma(fout, res, kappa, -1.5, "G_K_M3_2_H");

    mpfr_mul_ui(x, kappa, 2, RND);
    mpfr_gamma(res, x, RND);                   // x = gamma(2 * kappa)
    foutput(fout, res, "G_2K_H");

    calc_gamma_minus(fout, res, kappa, -0.5, "G_M1_2_MK_H");
    calc_gamma_minus(fout, res, kappa, 0.5, "G_1_2_MK_H");
    calc_gamma_minus(fout, res, kappa, 2.5, "G_5_2_MK_H");

    calc_lambda_vcj_p2(res, kappa, RHO_H, N0H_BY_N0E, x);
    foutput(fout, res, "LAMBDA_VC_P2_H");


    // Cold specie calculations
    const double rho_c = RHO_H / sqrt(TH_BY_TC);
    const double n0c_by_n0e = 1 - N0H_BY_N0E;

    fprintf(fout, "\n\n#define RHO_C %.17g", rho_c);            // .17g guarantees that the full double is printed
    fprintf(fout, "\n#define N0C_BY_N0E %.17g", n0c_by_n0e);

    mpfr_set_d(kappa, KAPPA_C, RND);

    calc_gamma(fout, res, kappa, 0.5, "G_K_P1_2_C");
    calc_gamma(fout, res, kappa, -0.5, "G_K_M1_2_C");
    calc_gamma(fout, res, kappa, 1, "G_K_P1_C");
    calc_gamma(fout, res, kappa, -1, "G_K_M1_C");
    calc_gamma(fout, res, kappa, 1.5, "G_K_P3_2_C");
    calc_gamma(fout, res, kappa, -1.5, "G_K_M3_2_C");

    mpfr_mul_ui(x, kappa, 2, RND);
    mpfr_gamma(res, x, RND);                   // x = gamma(2 * kappa)
    foutput(fout, res, "G_2K_C");

    calc_gamma_minus(fout, res, kappa, -0.5, "G_M1_2_MK_C");
    calc_gamma_minus(fout, res, kappa, 0.5, "G_1_2_MK_C");
    calc_gamma_minus(fout, res, kappa, 2.5, "G_5_2_MK_C");

    calc_lambda_vcj_p2(res, kappa, rho_c, n0c_by_n0e, x);
    foutput(fout, res, "LAMBDA_VC_P2_C");


    fprintf(fout, "\n");
    fclose(fout);


    mpfr_clears(res, kappa, x, (mpfr_t *) 0);
}


void foutput(FILE * fout, mpfr_t x, char * name)
{
    fprintf(fout, "\n#define %s \"", name);
    mpfr_out_str(fout, 10, 0, x, RND);
    fprintf(fout, "\"");
}


void calc_lambda_vcj_p2(mpfr_t res, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x)
{
    if (mpfr_inf_p(kappa))
    {
            mpfr_set_inf(res, 1);       // It is known from taking the limit to infinity of the expression for lambda_vcj_p2 that the answer is infinity
            return;
    }

    mpfr_set_d(res, 3 * LAMBDA + 1, RND);        // r = 3 * LAMBDA + 1

    mpfr_add_ui(x, kappa, 1, RND);
    mpfr_mul(res, res, x, RND);                  // r *= (kappa + 1)

    mpfr_add_d(x, kappa, 1.5, RND);
    mpfr_gamma(x, x, RND);
    mpfr_mul(res, res, x, RND);                 // r *= gamma(kappa + 1.5)

    mpfr_sub_d(x, kappa, 1.5, RND);
    mpfr_mul(res, res, x, RND);                 // r *= (kappa - 1.5)
    mpfr_mul_d(res, res, pow(rho, 2), RND);     // r *= pow(rho_j, 2)
    mpfr_div_d(res, res, n0_by_n0e, RND);       // r /= n0_by_n0e
    mpfr_div_d(res, res, pow(OMEGA_UH_BY_OMEGA_CE, 2) - 1, RND);    // r /= (pow(OMEGA_UH_BY_OMEGA_CE, 2)  - 1)

    mpfr_sub_d(x, kappa, 0.5, RND);
    mpfr_div(res, res, x, RND);                 // r /= (kappa - 0.5)
}


void calc_gamma(FILE * fout, mpfr_t res, mpfr_t kappa, double delta, char * name)
{
    if (delta > 0)
    {
        mpfr_add_d(res, kappa, delta, RND);
    }
    else
    {
        delta *= -1;
        mpfr_sub_d(res, kappa, delta, RND);
    }

    mpfr_gamma(res, res, RND);

    foutput(fout, res, name);
}


void calc_gamma_minus(FILE * fout, mpfr_t res, mpfr_t kappa, double delta, char * name)
{
    mpfr_set_d(res, delta, RND);            // res = delta
    mpfr_sub(res, res, kappa, RND);         // res -= kappa
    mpfr_gamma(res, res, RND);              // res = gamma(delta - kappa)

    foutput(fout, res, name);
}
