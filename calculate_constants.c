/*
 * Calculate the constants derived from the ones defined in constants.h and store
 * them in derived.h.
 *
 * This process will be carried out before compilation so that these values are
 * pre-compiled in to the executable and the calculations won't have to be
 * carried out over and over again.
 */

#include <stdio.h>
#include <math.h>
#include "constants.h"
#include <mpfr.h>


void calc_lambda_vcj_p2(mpfr_t result, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x);


int main(void)
{
    mpfr_set_default_prec(512);

    mpfr_t res, kappa, x;
    mpfr_inits(res, kappa, x, (mpfr_t *) 0);


    const double rho_c = RHO_H / sqrt(TH_BY_TC);
    const double n0c_by_n0e = 1 - N0H_BY_N0E;

    FILE * fout = fopen("derived.h", "w");

    fprintf(fout, "#define RHO_C %.17g", rho_c);            // .17g guarantees that the full double is printed
    fprintf(fout, "\n#define N0C_BY_N0E %.17g", n0c_by_n0e);


    // Hot specie calculations
    mpfr_set_d(kappa, KAPPA_H, RND);

    calc_lambda_vcj_p2(res, kappa, RHO_H, N0H_BY_N0E, x);
    fprintf(fout, "\n#define LAMBDA_VC_P2_H \"");
    mpfr_out_str(fout, 10, 0, res, RND);
    fprintf(fout, "\"");


    // Cold specie calculations
    mpfr_set_d(kappa, KAPPA_C, RND);

    calc_lambda_vcj_p2(res, kappa, rho_c, n0c_by_n0e, x);
    fprintf(fout, "\n#define LAMBDA_VC_P2_C \"");
    mpfr_out_str(fout, 10, 0, res, RND);
    fprintf(fout, "\"");


    fprintf(fout, "\n");
    fclose(fout);


    mpfr_clears(res, kappa, x, (mpfr_t *) 0);
}


void calc_lambda_vcj_p2(mpfr_t res, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x)
{
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
    mpfr_div(res, res, x, RND);                 // r =/ (kappa - 0.5)
}
