/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include "functions.h"


// Define a structure to store the start and end of the summation to calculate
// the dispersion value
struct Limits {
        int start;
        int end;
};

// function prototype of 'limits' function that calculates the summation limits
// for calculating the Dispersion relation
struct Limits limits(double omega);


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
 * Compute this part of the function separately since it has complicated
 * behaviour at x = 0
 */
double Gamma_n_by_x2(int n, double x)
{
        /*
         * From the definition of I_n(x) (http://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html, in particular the summation form (3)) we can prove that Gamma_n_by_x2 for n = 1 and x = 0 is 0.5 otherwise for n > 1 and x = 0 it is zero.
         */
        if (x == 0)
        {
                if (n == 1)
                        return 0.5;
                else
                        return 0.0;
        }

        return Gamma_n(n, x) / (x * x);
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
 * Define the summand inside the Dispersion relation which is function of the
 * index n, k_perp, and Omega.
 */
double Summand_n(int n, double k_perp, double omega)
{
        double beta_c = k_perp * RHO_C;
        double single = omega / (n * OMEGA_C);
        double denom = single * single - 1;

        return 2 * OMEGA_P_2 * Gamma_n(n, beta_c) / (beta_c * beta_c * OMEGA_C_2 * denom);
}


/*
 * Use the omega value to calculate the relevant limits for summation for the
 * dispersion relation.
 * Given the nature of the relation we want to sum the HALF_MAX_N terms on
 * either side of the value omega.
 */
struct Limits limits(double omega)
{
        struct Limits ans;

        ans.start = (int) omega - HALF_MAX_N;
        ans.end = (int) omega + HALF_MAX_N;

        // If the start is less than 1 we have to restrict it to 1 (the lowest)
        // index of the summation
        if (ans.start < 1)
                ans.start = 1;

        return ans;
}


/*
 * Define the disperstion relation as a function of k_perp and omega.
 */
double D(double k_perp, double omega)
{
        int n;
        double sum = 0;

        // Calculate the summation limits based on the value of omega
        struct Limits l = limits(omega);

        for (n = l.start; n <= l.end; n++)
                sum += Summand_n(n, k_perp, omega);

        return 1 - sum;
}


/*
 * Define function that calculates D() over an array of values.
 */
void D_array(double omega[], double Ds[], int size)
{
        int i;

        for (i = 0; i < size; i++)
                Ds[i] = D(K_PERP, omega[i]);

        return;
}


/*
 * Define a wrapper for D( , ) to make it more palatable for the gsl root finding
 * apparatus.
 */
// struct D_params {
//         double omega;
// };


double D_root(double k_perp, void *params)
{
        struct D_params *d_params = (struct D_params *) params;

        return D(k_perp, d_params->omega);
}
