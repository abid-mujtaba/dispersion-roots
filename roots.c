#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "dispersion.h"
#include "roots.h"


#define LOGFILEFMT "data/log-%d.txt"



// Declare the function D to refer to the Henning Dispersion relation/function
// double (* D_function)(double, double) = D_Henning;
double (* D_function)(double, double) = D;


/*
 * Use the bisection method for root finding to calculate the omega root of D( , ) for the specified
 * value of k_perp in the specified interval (lo to hi).
 *
 * If the root is found return 1 (True)
 * The found root is stored using the pointer 'root'
 */
int find_omega_root(const double k_perp, const double abs_lo, const double abs_hi, double * root, const double guess)
{
        double lo, hi, r, D_r;

        // Place the interval about the guessed value
        lo = guess - 1e-2;
        hi = guess + 1e-2;

        if (lo < abs_lo)
            lo = abs_lo;

        if (hi > abs_hi)
            hi = abs_hi;

        double D_lo = D_function(k_perp, lo);
        double D_hi = D_function(k_perp, hi);

        // If the guessed values don't span the root we revert to the absolute
        // values
        if (D_lo * D_hi > 0)
        {
            D_lo = D_function(k_perp, abs_lo);
            D_hi = D_function(k_perp, abs_hi);

            // We start by checking if the function changes signs at the absolute end-points. If they do we return 0 to indicate that a root could not be found
            if (D_lo * D_hi > 0)
                return 0;

            lo = abs_lo;
            hi = abs_hi;

        }


        while (hi - lo > ROOT_INTERVAL)
        {
            r = (lo + hi) / 2;
            D_r = D_function(k_perp, r);

            if (D_r * D_lo > 0)         // D has the same sign at r and lo so by bisection move lo to r
                lo = r;
            else
                hi = r;
        }

        *root = r;              // Store the result using the pointer 'root'

        return 1;
}


// Look for roots of D(,) at the k_perp sample values specified.

// Return value is the number of actual roots found (some might be out of the interval range)
int find_omega_roots_array(const int initial, double k_perp_samples[], double k_perps[], double omegas[], const int size)
{
        // Open log file for debug messages
        char logfile[80];       // Stores the constructed name of the log-file.
        FILE * fout = NULL;

        if (DEBUG)
        {
            sprintf(logfile, LOGFILEFMT, initial);
            fout = fopen(logfile, "w");
        }


        int count = 0;

        // Define the root finding intervals shifted from the ends (since the function D(,) is not defined there)
        double lo = initial + DELTA;
        double hi = initial + 1 - DELTA;
        double guess = initial + 0.5;    // First guess is the mid-point

        for (int i = 0; i < size; ++i)
        {
                k_perps[count] = k_perp_samples[i];

                if (DEBUG)
                    if (fout) 
                    {
                        fprintf(fout, "Finding root at k_perp = %.2f\n", k_perps[count]);
                        fflush(fout);
                    }

                if (find_omega_root(k_perp_samples[i], lo, hi, &omegas[count], guess))
                {
                        guess = omegas[count++];        // The guess is updated to be the last found root. It is hoped that continuity means the next root is close by
                }
        }

        if (fout)
            fclose(fout);

        return count;
}
