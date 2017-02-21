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
void I_n_array(int n, double x[], double In[], int size)
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


/*
 * Define the summand inside the Dispersion relation which is a function of the
 * integer n and Omega
 */
double Summand_n(int n, double omega)
{
 double single = omega / (n * OMEGA_C);
 double denom = single * single - 1;

 return 2 * OMEGA_P_2 * Gamma_n(n, BETA_C) / (BETA_C_2 * OMEGA_C_2 * denom);
}


/*
* Define the dispersion relation as a function of omega. The value of k_perp
* is fixed (in functions.h)
*/
double D(double omega)
{
        int n;
        double sum = 0;

        for (n = 1; n <= MAX_N; n++)
                sum += Summand_n(n, omega);

        return 1 - sum;
}


/*
 * Define function that calculates D() over an array of values.
 */
void D_array(double x[], double Ds[], int size)
{
        int i;

        for (i = 0; i < size; i++)
                Ds[i] = D(x[i]);

        return;
}
