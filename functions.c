/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include "functions.h"

/*
 * The second array is populated with the result of I_n (modified bessel fn)
 * applied to the values of the first array.
 */
void In_array(int n, double x[], double In[], int size)
{
        int i;

        for (i = 0; i < size; i++)
                In[i] = gsl_sf_bessel_In(n, x[i]);

        return;
}


/*
 * \Gamma_n(x) = \exp^{-x^2} I_n(x^2)
 * where I_n is the modified bessel function.
 */
double Gamma_n(int n, double x)
{
        double x2 = x * x;

        return gsl_sf_exp(-x2) * gsl_sf_bessel_In(n, x2);
}


/*
 * The second array is populated with the result of Gamma_n applied to the
 * values of the first array.
 */
void Gamma_n_array(int n, double x[], double Gn[], int size)
{
        int i;

        for (i = 0; i < size; i++)
                Gn[i] = Gamma_n(n, x[i]);

        return;
}
