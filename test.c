#include <stdio.h>
#include <math.h>
#include <mpfr.h>
#include "constants.h"
#include "functions.h"
#include "hypergeom.h"


int main(void)
{
        double k = 0;

        while (k <= 100)
        {
            double k_perp = k;
            double omega = 1.85;

            k += 10;

            // double two_lambda_H_prime = calc_two_lambda_j_prime(KAPPA_H, RHO_H, k_perp);

            /*mpfr_set_default_prec(PRECISION);*/

            /*mpfr_t om_by_om_cj, two_lambda_H_prime, coeff, h1f2, h2f3;*/
            /*mpfr_inits(om_by_om_cj, two_lambda_H_prime, coeff, h1f2, h2f3, (mpfr_ptr) 0);*/

            /*calc_omega_by_omega_cj(om_by_om_cj, omega, OMEGA_CH);*/
            /*calc_two_lambda_j_prime(two_lambda_H_prime, KAPPA_H, RHO_H, k_perp);*/


            /*struct coeffs_1f2 c_1f2;*/
            /*struct coeffs_2f3 c_2f3;*/

            /*init_coeffs(& c_1f2, & c_2f3);*/

            /*calc_coeffs_1f2(& c_1f2, KAPPA_H, om_by_om_cj);*/
            /*calc_coeffs_2f3(& c_2f3, KAPPA_H, om_by_om_cj);*/


            /*calc_coeff(coeff, om_by_om_cj, KAPPA_H, two_lambda_H_prime);*/
            /*hyp1F2(h1f2, c_1f2, two_lambda_H_prime);*/
            /*hyp2F3(h2f3, c_2f3, two_lambda_H_prime);*/

            printf("\n\nk_perp = %f", k_perp);
            /*mpfr_printf("\nomega_by_omega_cj = %RG", om_by_om_cj);*/
            /*mpfr_printf("\ntwo_lambda_j_prime = %RG", two_lambda_H_prime);*/
            /*mpfr_printf("\ncoeff = %RG", coeff);*/
            /*mpfr_printf("\n1F2(k_perp = %.1f, omega = %.1f) = %.RG", k_perp, omega, h1f2);*/
            /*mpfr_printf("\n2F3(k_perp = %.1f, omega = %.1f) = %.RG", k_perp, omega, h2f3);*/
            printf("\nD(%.1f, %.1f) = %f", k_perp, omega, D(k_perp, omega));

            /*mpfr_clears(om_by_om_cj, two_lambda_H_prime, coeff, h1f2, h2f3, (mpfr_ptr) 0);*/
            /*clear_coeffs(& c_1f2, & c_2f3);*/
        }


        printf("\n\n");

        return 0;
}
