#include <stdio.h>
#include "dispersion.h"
#include "constants.h"
#include <gsl/gsl_sf_gamma.h>


int main(void)
{
        double k_perp = 0.0;
        double omega = 1.6;

        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        printf("\n\n");

        return 0;
}
