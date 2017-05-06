#include <stdio.h>
#include "dispersion.h"
#include "constants.h"
#include <gsl/gsl_sf_gamma.h>


double lambda_vcj_p2(double, double, double);


int main(void)
{
        // double k_perp = 0.1;
        // double omega = 0.6;
        //
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));
        //
        // omega = 1.6;
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));
        //
        // omega = 2.5;
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));
        //
        // k_perp = 1.4;
        // omega = 3.6;
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));
        //
        // omega = 4.6;
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));
        //
        // omega = 5.5;
        // printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));


        double l = lambda_vcj_p2(KAPPA_H, RHO_H, N0H_BY_N0E);
        double denom = (KAPPA_H - 1.5) * (3 * LAMBDA + 1) * (KAPPA_H + 1) * (KAPPA_H + 0.5) * gsl_sf_gamma(KAPPA_H - 0.5);
        double b2 = denom / l;

        printf("\nlambda_vcj_p2 = %.17g", l);
        printf("\ndenom = %.17g", denom);
        printf("\nbeta^2 = %.17g", b2);


        printf("\n\n");

        return 0;
}
