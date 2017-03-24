#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "functions.h"
#include "hypergeom.h"


int main(void)
{
        // int i;
        //
        // for (i = 0; i < 100; ++i)
        // {
        //         double k = i / 5.0;
        //
        //         printf("\nD(%.2f, 0.5) = %.5e", k, D(k, 0.5));
        // }
        //
        // printf("\n\n");
        //
        // for (i = 0; i < 20; ++i)
        // {
        //         double w = i / 20.0;
        //
        //         printf("\nD(15, %.2f) = %.5e", w, D(15, w));
        // }

        // double k = 6.8, w = 0.5;
        // printf("\n\nD(%.1f, %.1f) = %.20f", k, w, D(k, w));

        double k_perp = 6.8;
        double omega = 0.5;

        double omega_by_omega_cj = omega / OMEGA_CH;
        double two_lambda_j_prime = calc_two_lambda_j_prime(KAPPA_H, RHO_H, k_perp);
        struct coeffs_1f2 c_1f2 = calc_coeffs_1f2(KAPPA_H, omega_by_omega_cj);

        printf("\n1F2(k_perp = %.1f, omega = %.1f) = %.17g", k_perp, omega, hyp1F2(c_1f2, two_lambda_j_prime));

        printf("\n\n");

        return 0;
}
