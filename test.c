#include <stdio.h>
#include "dispersion.h"
#include "constants.h"
#include <gsl/gsl_sf_gamma.h>


double lambda_vcj_p2(double, double, double);


int main(void)
{
        double k_perp = 0.1;
        double omega = 0.6;

        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        omega = 1.6;
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        omega = 2.5;
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        k_perp = 1.4;
        omega = 3.6;
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        omega = 4.6;
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        omega = 5.5;
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));


        printf("\n\n");

        return 0;
}
