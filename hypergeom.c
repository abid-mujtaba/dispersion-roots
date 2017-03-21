/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hypergeom.h"


double hyp1F2(const struct coeffs_1f2 c, const double x)
{
        if (x == 0)                     // Base case of recursion
                return 1;


        int k, i;

        double result = 0;
        double terms[MAX_TERMS];
        double term = 1;                // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE

        // We split the value of x in to a number of steps and the remainder for incremental calculation using power Series
        double remainder = fmod(x, TAYLOR_STEP);
        double step = x - remainder;

        // If we are exactly on a step we must move to the previous step
        if (remainder == 0)
        {
                remainder = TAYLOR_STEP;
                step -= TAYLOR_STEP;
        }

        // Calculate first term in Taylor Series about x = THRESHOLD
        double coeff = 1;
        struct coeffs_1f2 c2;

        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                c2.a1 = c.a1 + k;
                c2.b1 = c.b1 + k;
                c2.b2 = c.b2 + k;

                terms[k] = coeff * hyp1F2(c2, step);

                coeff *= c2.a1 * remainder / (c2.b1 * c2.b2 * (k + 1));
        }

        printf("\rx = %f - k = %d", x, k);

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

                terms[k] = coeff * hyp2F3(c2, step);

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
