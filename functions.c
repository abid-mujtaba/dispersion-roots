/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <gsl/gsl_sf_bessel.h>
#include "functions.h"

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
