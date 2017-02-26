/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>              // Defines GSL_SUCCESS and GSL_CONTINUE
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

        return 2 * OMEGA_P_2 * Gamma_n_by_x2(n, beta_c) / (OMEGA_C_2 * denom);
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
 * apparatus. This also requires a struct to carry additional parameters for
 * the function, in this case the value of omega since we want to find the value
 * of k_perp that is a root of D( , ) for a specified value of omega.
 */
struct D_params {
        double omega;           // the only additional value needed (apart for k_perp) to calculate D( , )
};


double D_root(double k_perp, void *params)
{
        // Declare a struct D_params pointer object and use type-casting to convert void *params to this type
        struct D_params *d_params = (struct D_params *) params;

        // Extract the value of omega from the pointer (using ->) and send it to the D( , ) function
        return D(k_perp, d_params->omega);
}


/*
 * Use gsl root finding to calculate the k_perp root of D( , ) for the specified
 * value of omega in the specified interval (lo to hi).
 */
double find_k_perp_root(double omega, double lo, double hi)
{
        // Create D_params and store specified value of omega inside it
        struct D_params params;
        params.omega = omega;

        // Create a gsl_function object and store the function whose root needs to be fount (D_root which has the correct signature, note the void *params) and the params that said function needs
        gsl_function F;
        F.function = &D_root;
        F.params = &params;

        // Create and populate the solver and solver_type
        const gsl_root_fsolver_type *solver_type;
        gsl_root_fsolver *solver;

        solver_type = gsl_root_fsolver_brent;                   // Using the Brent method
        solver = gsl_root_fsolver_alloc(solver_type);           // Create the solver. This alloc needs a corresponding 'free' call at end

        // Tell gsl_root about the solver to use, the function whose root is to be found and the bracket limits for finding the root
        gsl_root_fsolver_set(solver, &F, lo, hi);


        double r, low, high;
        int i, test_status = GSL_CONTINUE;

        for (i = 0; i <= ROOT_MAX_ITERATIONS && test_status == GSL_CONTINUE; ++i)
        {
                // Iterate the solver once. This will check the function at the bracket limits, calculate a new estimate of the root and narrow down the bracket
                // It returns a status/error code to indicate how the iteration went. It doesn't comment on whether convergence was achieved
                gsl_root_fsolver_iterate(solver);

                // Get values of interest after the iteration
                r = gsl_root_fsolver_root(solver);
                low = gsl_root_fsolver_x_lower(solver);
                high = gsl_root_fsolver_x_upper(solver);

                // We test the situation after the iteration by comparing the new bracket with the required ROOT_INTERVAL
                // The return value is a status/error code which will indicate success or the need to continue
                test_status = gsl_root_test_interval(low, high, 0, ROOT_INTERVAL);

                if (test_status == GSL_SUCCESS)                 // Root has been found to within specified interval
                        break;
        }

        // Free the solver once we no longer need it
        gsl_root_fsolver_free(solver);

        if (test_status != GSL_SUCCESS)
                printf("\nUnable to find root for omega = %.2f. Root test status: %d = %s.", omega, test_status, gsl_strerror(test_status));

        return r;
}

/*
 * Take an array of omega values (slices) where the corresponding k_perp root is
 * to be found. Since not every omega has a root we will returned modified arrays
 * containing the omega and k_perp values where a root was found.
 *
 * The function will return an int representing the size of the returned arrays.
 */
int find_k_perp_roots_array(double slices[], double omega[], double roots[], int size)
{
        int i, j, count = 0;
        double om;
        double mid = (ROOT_LO + ROOT_HI) * 0.5;               // Midway point for dealing with two roots
        double dLo, dHi;

        for (i = 0; i < size; ++i)
        {
                om = slices[i];

                /* first we check if their is a sign-change for the bracket limits.
                 * if this is not the case then the root multiplicity is even which might mean zero or two roots (for this function)
                 * Note: signbit returns non-zero (128) if its argument's sign bit is set (negative number)
                 *
                 * We use the bit-wise XOR operator (^) to determine if the function value at the bracket limits has different signs which is neccessary for the root-finding to work
                 */

                dLo = D(ROOT_LO, om);
                dHi = D(ROOT_HI, om);

                if (signbit(dLo) ^ signbit(dHi))
                {
                        omega[count] = om;
                        roots[count++] = find_k_perp_root(om, ROOT_LO, ROOT_HI);
                }
                else            // Zero or two roots
                {
                        /*
                         * If ROOT_LO and ROOT_HI do not have opposite polarity then we are possibly dealing with a case of two roots.
                         * We set up a for loop to look for a possible mid value that splits the range in to two containing roots in each.
                         * The for loop will iterate up to a maximum of 4 times.
                         * Note the test in the middle which tests both j < 4 (to limit the iterations) and the polarity at lo and mid.
                         * If at the beginning the previous (heuristic value of mid) has a flip the loop will NOT even be initiated
                         * If it is initiated then the first value will be midway between LO and HI
                         */
                        for (j = 0; j < 4 && !(signbit(dLo) ^ signbit(D(mid, om))); ++j)
                                mid = ROOT_LO + (ROOT_HI - ROOT_LO) / pow(2, j + 1);            // Previous value of mid didn't work so move it closer to ROOT_LO in a binary fashion.

                        if (signbit(dLo) ^ signbit(D(mid, om)))         // Flip mid-way so two roots found
                        {
                                omega[count] = om;
                                roots[count++] = find_k_perp_root(om, ROOT_LO, mid);

                                omega[count] = om;
                                roots[count++] = find_k_perp_root(om, mid, ROOT_HI);

                                /*
                                 * We are aware that the roots at the slice above the current one will be close to these roots so we set mid to be mid-way between these roots to speed up the next convergence.
                                 */
                                mid = (roots[count - 2] + roots[count - 1]) * 0.5;
                        }
                        else {                                          // Zero roots. Reset value of mid
                                mid = (ROOT_LO + ROOT_HI) * 0.5;
                        }
                }
        }

        return count;
}
