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

// Take two arrays with specified length
// Triple every value in the first array and store it in the second array to
// return it by reference
void triple(double a[], double t[], int n)
{
        int i;

        for (i = 0; i < n; i++)
                t[i] = a[i] * 3;

        return;
}


/*
 * Take two arrays. The first is a list of x-values.
 * The second will carry the result back by reference.
 * The result is J0(.) applied to every value of the first.
 */
void J0_array(double x[], double j0[], int size)
{
        int i;

        for (i = 0; i < size; i++)
                j0[i] = gsl_sf_bessel_J0(x[i]);

        return ;
}
