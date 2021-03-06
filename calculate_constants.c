/*
 * Calculate the constants derived from the ones defined in constants.h and store
 * them in derived.h.
 *
 * This process will be carried out before compilation so that these values are
 * pre-compiled in to the executable and the calculations won't have to be
 * carried out over and over again.
 */

 //       Use recursion where necessary to reduce calculations

#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "roots.h"
#include <mpfr.h>


// Store the mpfr var x with name in the output file:  #define <name> "<value of x>"
void foutput(FILE * fout, mpfr_t x, char * name);

void calc_lambda_vcj_p2(mpfr_t result, mpfr_t lambda, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x, mpfr_t y);


int main(void)
{
    // Calculate the max needed precision from max value of k_perp
    // All constants are calculated with this max precision since downgrading is easy
    int p = 1 + (int) (K_PERP_MAX / 30);
    mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

    mpfr_t res, kappa, lambda, x, y;
    mpfr_inits(res, kappa, lambda, x, y, (mpfr_t *) 0);

    FILE * fout = fopen("derived.h", "w");


    // Hot specie calculations
    mpfr_set_d(kappa, KAPPA_H, RND);
    mpfr_set_d(lambda, LAMBDA_H, RND);

    calc_lambda_vcj_p2(res, lambda, kappa, RHO_H, N0H_BY_N0E, x, y);
    foutput(fout, res, "LAMBDA_VC_P2_H");


    // Cold specie calculations
    const double rho_c = RHO_H / sqrt(TH_BY_TC);
    const double n0c_by_n0e = 1 - N0H_BY_N0E;

    fprintf(fout, "\n\n#define RHO_C %.17g", rho_c);            // .17g guarantees that the full double is printed
    fprintf(fout, "\n#define N0C_BY_N0E %.17g", n0c_by_n0e);

    mpfr_set_d(kappa, KAPPA_C, RND);
    mpfr_set_d(lambda, LAMBDA_C, RND);

    calc_lambda_vcj_p2(res, lambda, kappa, rho_c, n0c_by_n0e, x, y);
    foutput(fout, res, "LAMBDA_VC_P2_C");


    fprintf(fout, "\n");
    fclose(fout);


    mpfr_clears(res, kappa, lambda, x, (mpfr_t *) 0);
}


void foutput(FILE * fout, mpfr_t x, char * name)
{
    fprintf(fout, "\n#define %s \"", name);
    mpfr_out_str(fout, 10, 0, x, RND);
    fprintf(fout, "\"");
}


// lambda^2_vcj = v^2_th / w^2_pj (1 + 3 Lambda (k - 3/2) / (k - 5/2)) (k - 3/2) (k + 1/2) (k + 1)
// where v^2_th / w^2_pj = rho^2_j / (n0j / n0e * ((w_UH / w_ce)^2 - 1))
//
// For the case kappa -> infinity all the kappa terms become 1 since the kappa^3 order of lambda_vcj_p2 is equal to the kappa^3 order of the alpha[n]
void calc_lambda_vcj_p2(mpfr_t res, mpfr_t lambda, mpfr_t kappa, double rho, double n0_by_n0e, mpfr_t x, mpfr_t y)
{
    mpfr_mul_ui(y, lambda, 15, RND);				// y = lambda * 15
    mpfr_div_ui(y, y, 4, RND);						// y = lambda * 15 / 4
    mpfr_set_ui(res, 1, RND);

    if (! mpfr_inf_p(kappa))
    {
        mpfr_sub_d(x, kappa, 1.5, RND);
        mpfr_mul(res, res, x, RND);                 // r *= (k - 3/2)

        mpfr_mul(y, y, x, RND);						// y = 15 * Lambda / 4 * (k - 3/2)

        mpfr_sub_d(x, kappa, 2.5, RND);
        mpfr_div(y, y, x, RND);						// y = 15 * Lambda / 4 * (k - 3/2) / (k - 5/2)

        mpfr_add_d(x, kappa, 0.5, RND);
        mpfr_mul(res, res, x, RND);                 // r *= (k + 1/2)

        mpfr_add_ui(x, kappa, 1, RND);
        mpfr_mul(res, res, x, RND);                 // r *= (k + 1)
    }

	mpfr_add_ui(y, y, 1, RND);					// y = 1 + 15 / 4 * Lambda * (k - 3/2) / (k - 5/2)	(if k != inf)
    mpfr_mul(res, res, y, RND);					// r *= y 

    mpfr_mul_d(res, res, pow(rho, 2), RND);     // r *= pow(rho_j, 2)

    mpfr_div_d(res, res, n0_by_n0e, RND);       // r /= n0_by_n0e
    mpfr_div_d(res, res, pow(OMEGA_UH_BY_OMEGA_CE, 2) - 1, RND);    // r /= (pow(OMEGA_UH_BY_OMEGA_CE, 2)  - 1)
}
