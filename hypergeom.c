/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hypergeom.h"


double hyp1F2(const struct coeffs_1f2 c, const double x, double steps_1F2[NUM_STEPS][MAX_TERMS], int populating)
{
        if (x == 0)                     // Base case of recursion
                return 1;


        int k, i;

        double result = 0;
        double terms[MAX_TERMS];
        double term = 1;                // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE

        // We split the value of x in to a number of steps and the remainder for incremental calculation using power Series
        double remainder = fmod(x, TAYLOR_STEP);
        int step = (x - remainder) / TAYLOR_STEP;

        // If we are exactly on a step, we already have the answer so we can return it directly
        if (remainder == 0)
        {
                if (! populating)
                        return steps_1F2[step][0];              // All higher terms will go to zero when remainder = 0 (in the expansion)

                // If 'populating' is True we need to calculate the values for which we move down one step perform the calculation using remainder = TAYLOR_STEP

                step -= 1;
                remainder = TAYLOR_STEP;
        }

        double coeff = 1;
        struct coeffs_1f2 c2 = c;

        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                term = coeff * steps_1F2[step][k];
                terms[k] = term;

                coeff *= c2.a1 * remainder / (c2.b1 * c2.b2 * (k + 1));

                c2.a1 += 1;
                c2.b1 += 1;
                c2.b2 += 1;
        }

        // Sort the terms before adding them to reduce computational errors from adding disparate numbers
        qsort(terms, k, sizeof(double), compare_terms);

        for (i = 0; i < k; ++i)
                result += terms[i];

        return result;
}


double hyp2F3(const struct coeffs_2f3 c, const double x)
{
        if (x == 0)             // Base case of recursion
                return 1;


        int k, i;

        double result = 0;
        double terms[MAX_TERMS];
        double term = 1;


        double remainder = fmod(x, TAYLOR_STEP);
        double step = x - remainder;

        if (remainder == 0)
        {
                remainder = TAYLOR_STEP;
                step -= TAYLOR_STEP;
        }


        double coeff = 1;
        struct coeffs_2f3 c2;

        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                c2.a1 = c.a1 + k;
                c2.a2 = c.a2 + k;
                c2.b1 = c.b1 + k;
                c2.b2 = c.b2 + k;
                c2.b3 = c.b3 + k;

                term = coeff * hyp2F3(c2, step);
                terms[k] = term;

                coeff *= c2.a1 * c2.a2 * remainder / (c2.b1 * c2.b2 * c2.b3 * (k + 1));
        }


        qsort(terms, k, sizeof(double), compare_terms);

        for (i = 0; i < k; ++i)
                result += terms[i];

        return result;
}


/*
 * Define a function for comparing terms in the array based on their absolute value.
 */
int compare_terms(const void *pa, const void *pb)
{
        // Use the void * pointers to get the long double values and take their absolute value.
        double a = fabsl(* (double *) pa);
        double b = fabsl(* (double *) pb);

        // The return values are chosen to give us ascending order
        if (a < b)
                return -1;

        if (a > b)
                return 1;

        return 0;
}


/*
 * Populate the intermediate step values for 1F2 in the Taylor Expansion.
 */
void populate_steps_1F2(const struct coeffs_1f2 c, double steps_1F2[NUM_STEPS][MAX_TERMS])
{
        // Start by populating the zeroeth step where 1F2 = 0 regardless of the coefficients
        int i, k;

        for (i = 0; i < NUM_STEPS; ++i)
        {
                struct coeffs_1f2 c2 = c;

                for (k = 0; k < MAX_TERMS; ++k)
                {
                        steps_1F2[i][k] = hyp1F2(c2, TAYLOR_STEP * i, steps_1F2, 1);

                        c2.a1 += 1;
                        c2.b1 += 1;
                        c2.b2 += 1;
                }
        }
}
