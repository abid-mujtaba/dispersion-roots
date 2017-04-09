#include <stdio.h>
#include <math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>              // Defines GSL_SUCCESS and GSL_CONTINUE
#include "functions.h"
#include "roots.h"

/*
 * Define a wrapper for D( , ) to make it more palatable for the gsl root finding
 * apparatus. This also requires a struct to carry additional parameters for
 * the function, in this case the value of omega since we want to find the value
 * of k_perp that is a root of D( , ) for a specified value of omega.
 */
struct D_params {
        double omega;           // the only additional value needed (apart for k_perp) to calculate D( , )
        double k_perp;          // like-wise for calculating omega root. Only one is used at a time
};


double D_k_perp_root(const double k_perp, void *params)
{
        // Declare a struct D_params pointer object and use type-casting to convert void *params to this type
        struct D_params *d_params = (struct D_params *) params;

        // Extract the value of omega from the pointer (using ->) and send it to the D( , ) function
        return D(k_perp, d_params->omega);
}


double D_omega_root(const double omega, void *params)
{
        // Declare a struct D_params pointer object and use type-casting to convert void *params to this type
        struct D_params *d_params = (struct D_params *) params;

        // Extract the value of omega from the pointer (using ->) and send it to the D( , ) function
        return D(d_params->k_perp, omega);
}


/*
 * Use gsl root finding to calculate the k_perp root of D( , ) for the specified
 * value of omega in the specified interval (lo to hi).
 */
double find_k_perp_root(const double omega, const double lo, const double hi)
{
        // Create D_params and store specified value of omega inside it
        struct D_params params;
        params.omega = omega;

        // Create a gsl_function object and store the function whose root needs to be fount (D_root which has the correct signature, note the void *params) and the params that said function needs
        gsl_function F;
        F.function = &D_k_perp_root;
        F.params = &params;

        // Create and populate the solver and solver_type
        const gsl_root_fsolver_type *solver_type;
        gsl_root_fsolver *solver;

        solver_type = gsl_root_fsolver_brent;                   // Using the Brent method
        solver = gsl_root_fsolver_alloc(solver_type);           // Create the solver. This alloc needs a corresponding 'free' call at end

        // Tell gsl_root about the solver to use, the function whose root is to be found and the bracket limits for finding the root
        gsl_root_fsolver_set(solver, &F, lo, hi);


        double r, low, high;
        int test_status = GSL_CONTINUE;

        for (int i = 0; i <= ROOT_MAX_ITERATIONS && test_status == GSL_CONTINUE; ++i)
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
int find_k_perp_roots_array(double slices[], double omega[], double roots[], const int size)
{
        int count = 0;
        double om;
        double mid = (ROOT_LO + ROOT_HI) * 0.5;               // Midway point for dealing with two roots
        double dLo, dHi;

        for (int i = 0; i < size; ++i)
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
                        for (int j = 0; j < 4 && !(signbit(dLo) ^ signbit(D(mid, om))); ++j)
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


/*
 * Use gsl root finding to calculate the omega root of D( , ) for the specified
 * value of k_perp in the specified interval (lo to hi).
 *
 * If the root is found return 1 (True)
 * The found root is stored using the poiner 'root'
 */
int find_omega_root(const double k_perp, const double lo, const double hi, double * root)
{
        // We start by checking if the function changes signs at the end-points
        if (D(k_perp, lo) * D(k_perp, hi) > 0)          // If the two end-points have the same sign then this is true
                return 0;

        // Create D_params and store specified value of k_perp inside it
        struct D_params params;
        params.k_perp = k_perp;

        // Create a gsl_function object and store the function whose root needs to be fount (D_root which has the correct signature, note the void *params) and the params that said function needs
        gsl_function F;
        F.function = &D_omega_root;
        F.params = &params;

        // Create and populate the solver and solver_type
        const gsl_root_fsolver_type *solver_type;
        gsl_root_fsolver *solver;

        solver_type = gsl_root_fsolver_brent;                   // Using the Brent method
        solver = gsl_root_fsolver_alloc(solver_type);           // Create the solver. This alloc needs a corresponding 'free' call at end

        // Tell gsl_root about the solver to use, the function whose root is to be found and the bracket limits for finding the root
        gsl_root_fsolver_set(solver, &F, lo, hi);


        double r, low, high;
        int test_status = GSL_CONTINUE;

        for (int i = 0; i <= ROOT_MAX_ITERATIONS && test_status == GSL_CONTINUE; ++i)
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
                printf("\nUnable to find root for k_perp = %.2f. Root test status: %d = %s.", k_perp, test_status, gsl_strerror(test_status));

        *root = r;              // Store the result using the pointer 'root'

        return 1;
}


// Return value is the number of actual roots found (some might be out of the interval range)
int find_omega_roots_array(const int initial, double k_perps[], double omegas[], const int size)
{
        int count = 0;

        // Define the root finding intervals shifted from the ends (since the function D(,) is not defined there)
        double lo = initial + 1e-10;
        double hi = initial + 1 - 1e-10;
        double delta = K_PERP_MAX * 1.0 / size;

        k_perps[count] = 1e-5;             // The function D(,) does not change sign on k_perp = 0 so we edge a little forward
        if (find_omega_root(k_perps[count], lo, hi, &omegas[count]))
                ++count;

        for (int i = 1; i < size; ++i)
        {
                k_perps[count] = delta * i;

                printf("\r%d < omega < %d - Searching for root at k_perp = %.2f", initial, initial + 1, k_perps[count]);
                fflush(stdout);

                if (find_omega_root(k_perps[count], lo, hi, &omegas[count]))
                        ++count;
                else                            // This means the graph has fallen below our threshold. No point in searching further.
                        break;
        }

        return count;
}
