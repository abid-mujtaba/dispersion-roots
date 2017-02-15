/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int main(void)
{
        double x = 5.0;
        double y = gsl_sf_bessel_J0(x);

        printf("J0(%g) = %.18e\n", x, y);

        return 0;
}
