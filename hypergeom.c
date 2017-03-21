/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hypergeom.h"


double hyp1F2(const struct coeffs_1f2 c, const double x)
{
        int k, i;

        double result = 0;
        double terms[MAX_TERMS];
        double term = 1;                // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE

        // For small values of x we calculate 1F2 directly using the defining sum
        if (x <= THRESHOLD)
        {
                // When (the absolute value of) 'term' becomes less than TOLERANCE the sum stops changing rapidly and we truncate it
                for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
                {
                        terms[k] = term;

                        // Based on the definition of 1F2 each term differs from the previous by multiplicative factors that have to do with the Pochhammer symbols, power of x and the factorial in the denominator
                        term *= (c.a1 + k) * x / ((c.b1 + k) * (c.b2 + k) * (k + 1));
                }
        }
        else            // For large values of x use a Taylor Expansion about x = 4 for better convergence
        {
                // Calculate first term in Taylor Series about x =4
                double coeff = 1;
                struct coeffs_1f2 c2;

                for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
                {
                        c2.a1 = c.a1 + k;
                        c2.b1 = c.b1 + k;
                        c2.b2 = c.b2 + k;

                        terms[k] = coeff * hyp1F2(c2, 4);

                        coeff *= c2.a1 * (x - THRESHOLD) / (c2.b1 * c2.b2 * (k + 1));
                }
        }

        // Sort the terms before adding them to reduce computational errors from adding disparate numbers
        qsort(terms, k, sizeof(double), compare_terms);

        for (i = 0; i < k; ++i)
                result += terms[i];

        return result;
}


double hyp2F3(const struct coeffs_2f3 c, const double x)
{
        int k, i;

        double result = 0;
        double terms[MAX_TERMS];
        double term = 1;


        if (x <= THRESHOLD)
        {
                for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
                {
                        terms[k] = term;

                        term *= (c.a1 + k) * (c.a2 + k) * x / ((c.b1 + k) * (c.b2 + k) * (c.b3 + k) * (k + 1));
                }
        }
        else
        {
                double coeff = 1;
                struct coeffs_2f3 c2;

                for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
                {
                        c2.a1 = c.a1 + k;
                        c2.a2 = c.a2 + k;
                        c2.b1 = c.b1 + k;
                        c2.b2 = c.b2 + k;
                        c2.b3 = c.b3 + k;

                        terms[k] = coeff * hyp2F3(c2, 4);

                        coeff *= c2.a1 * c2.a2 * (x - THRESHOLD) / (c2.b1 * c2.b2 * c2.b3 * (k + 1));
                }
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
