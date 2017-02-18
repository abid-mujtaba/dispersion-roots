/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <gsl/gsl_sf_bessel.h>
#include "functions.h"

// Simply return the value from the gsl definition of this function
double J0(double x)
{
        return gsl_sf_bessel_J0(x);
}


// Calculate the sum of an array of doubles
double sum(double a[], int n)
{
        double result = 0;
        int i;

        for (i = 0; i < n; i++)
                result += a[i];

        return result;
}
