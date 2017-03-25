#include <stdio.h>
#include <math.h>
#include <mpfr.h>
#include "constants.h"
#include "functions.h"
#include "hypergeom.h"


int main(void)
{
        // int i;

        // for (i = 0; i < 100; ++i)
        // {
        //         double k = i / 5.0;
        //
        //         printf("\nD(%.2f, 0.5) = %.5e", k, D(k, 0.5));
        // }

        // printf("\n\n");

        // double k = 7;
        //
        // for (i = 0; i < 20; ++i)
        // {
        //         double w = i / 20.0;
        //
        //         printf("\nD(%.2f, %.2f) = %.5e", k, w, D(k, w));
        // }


        double k = 7.0, w = 0.85;
        printf("\n\nD(%.1f, %.1f) = %.20f", k, w, D(k, w));


        // double k_perp = 6.8;
        // double omega = 0.5;
        //
        // double omega_by_omega_cj = omega / OMEGA_CH;
        // double two_lambda_H_prime = calc_two_lambda_j_prime(KAPPA_H, RHO_H, k_perp);
        // struct coeffs_1f2 c_1f2 = calc_coeffs_1f2(KAPPA_H, omega_by_omega_cj);
        // struct coeffs_2f3 c_2f3 = calc_coeffs_2f3(KAPPA_H, omega_by_omega_cj);
        //
        // mpfr_set_default_prec(PRECISION);
        //
        // mpfr_t coeff, h1f2, h2f3;
        // mpfr_inits(coeff, h1f2, h2f3, (mpfr_ptr) 0);
        //
        // calc_coeff(coeff, omega_by_omega_cj, KAPPA_H, two_lambda_H_prime);
        // hyp1F2(h1f2, c_1f2, two_lambda_H_prime);
        // hyp2F3(h2f3, c_2f3, two_lambda_H_prime);
        //
        // mpfr_printf("\n\ncoeff = %RG", coeff);
        // mpfr_printf("\n1F2(k_perp = %.1f, omega = %.1f) = %.RG", k_perp, omega, h1f2);
        // mpfr_printf("\n2F3(k_perp = %.1f, omega = %.1f) = %.RG", k_perp, omega, h2f3);
        // printf("\nD(%.1f, %.1f) = %f", k_perp, omega, D(k_perp, omega));
        //
        // mpfr_clears(coeff, h1f2, h2f3, (mpfr_ptr) 0);


        printf("\n\n");

        return 0;
}
