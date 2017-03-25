#include <stdio.h>
#include <math.h>
#include <mpfr.h>
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


        // double k_perp = 6.8;
        // double omega = 0.5;
        //
        // double omega_by_omega_cj = omega / OMEGA_CH;
        // double two_lambda_j_prime = calc_two_lambda_j_prime(KAPPA_H, RHO_H, k_perp);
        // struct coeffs_1f2 c_1f2 = calc_coeffs_1f2(KAPPA_H, omega_by_omega_cj);
        //
        // printf("\n1F2(k_perp = %.1f, omega = %.1f) = %.17g", k_perp, omega, hyp1F2(c_1f2, two_lambda_j_prime));


        mpfr_set_default_prec(PRECISION);               // Set default precision for all variables whose precision is NOT explicitly specified when initialized

        mpfr_t x, y;               // Create a MPFR (float) variables
        mpfr_inits(x, y, (mpfr_ptr) 0);           // Initialize the variables. Note that the list of variables must be terminated with a NULL pointer of the correct type. Precision is NOT specified so the defaut value is used.
        mpfr_set_d(x, 3.14, MPFR_RNDN);                 // Set the value of 'x' to be equal to the specified double and use 'Nearest' rounding
        mpfr_mul_d(y, x, (1.0 / 3), MPFR_RNDN);

        mpfr_printf("\nx = %RG", x);
        mpfr_printf("\ny = %Rg", y);             // Print the MPFR variable (requires mpfr_print and the Rf specified)

        mpfr_printf("\n\nPrecision of x = %Pu bits", mpfr_get_prec(x));

        mpfr_clears(x, y, (mpfr_ptr) 0);                  // Clear the variable

        printf("\n\n");

        return 0;
}
