/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <math.h>
#include "hypergeom.h"


long double next_term_1F2(const double coeff, const struct coeffs_1f2 c_1f2, double x);
long double next_term_2F3(const struct coeffs_2f3 c_2f3, double x);
int terms_2F3(const struct coeffs_2f3, const double x, long double terms[]);


long double hyp1F2(const double a1, const double b1, const double b2, const double x)
{
        int k;


        long double result = 0;
        long double term = 1;                // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE

        // When (the absolute value of) 'term' becomes less than TOLERANCE the sum stops changing rapidly and we truncate it
        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                result += term;

                // Based on the definition of 1F2 each term differs from the previous by multiplicative factors that have to do with the Pochhammer symbols, power of x and the factorial in the denominator
                term *= (a1 + k) * x / ((b1 + k) * (b2 + k) * (k + 1));
        }

        return result;
}


long double hyp2F3(const double a1, const double a2, const double b1, const double b2, const double b3, const double x)
{
        int k;

        long double result = 0;
        long double term = 1;

        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                result += term;

                term *= (a1 + k) * (a2 + k) * x / ((b1 + k) * (b2 + k) * (b3 + k) * (k + 1));
        }

        return result;
}


long double together_hyp(const double coeff, const struct coeffs_1f2 c_1f2, const struct coeffs_2f3 c_2f3, const double x)
{
        int k;

        long double t_2f3[2 * MAX_TERMS];

        int size = terms_2F3(c_2f3, x, t_2f3);

        long double result = coeff * hyp1F2(c_1f2.a1, c_1f2.b1, c_1f2.b2, x);

        for (k = 0; k < size; ++k)
                result += t_2f3[k];

        return result;
}


long double series_hyp(const double coeff, const struct coeffs_1f2 c_1f2, const struct coeffs_2f3 c_2f3, const double x)
{
        int k = 0;

        long double result = 0;
        long double term_1f2 = coeff;           // The 1F2 is multiplied by the coefficient so we mutiply it directly with the initial term for 1F2
        long double term_2f3 = 1;
        long double term = term_1f2 - term_2f3;

        // Implemented using a do-while because if coeff == 1 then term = 0 initially but we still want the loop to execute at least once
        do
        {
                result += term;
                // printf("\nk = %2d  -  term_1f2 = %+.4Le  -  term_2f3 = %+.4Le  -  frac = %+.4Le  -  result = %+.4Le", k, term_1f2, term_2f3, term_1f2 / term_2f3, result);

                term_1f2 = next_term_1F2(coeff, c_1f2, x);
                term_2f3 = next_term_2F3(c_2f3, x);

                term = term_1f2 - term_2f3;

                ++k;            // Only use is to determine how many iterations have occurred
        }
        while ((k < MAX_TERMS) & (fabs(term) > TOLERANCE));

        return result;
}


long double next_term_1F2(const double coeff, const struct coeffs_1f2 c_1f2, const double x)
{
        // Note the use of static variables to maintain a running value of the summation index k and the value of the term itself
        static int k = 0;
        static long double term;

        if (k == 0)
                term = coeff;                   // Initial value of term

        term *= (c_1f2.a1 + k) * x / ((c_1f2.b1 + k) * (c_1f2.b2 + k) * (k + 1));

        ++k;

        return term;
}


/*
 * Calculate all of the terms of -2F3 and store them in the specified array
 */
int terms_2F3(const struct coeffs_2f3 c, const double x, long double terms[])
{
        int k;

        long double term = -1;          // - 2F3 is to be added to the final answer so each term will be negative

        for (k = 0; ((k < MAX_TERMS) & (fabs(term) > TOLERANCE)); ++k)
        {
                terms[k] = term;

                term *= (c.a1 + k) * (c.a2 + k) * x / ((c.b1 + k) * (c.b2 + k) * (c.b3 + k) * (k + 1));
        }

        return k;
}


long double next_term_2F3(const struct coeffs_2f3 c_2f3, const double x)
{
        static int k = 0;
        static long double term = 1;

        term *= (c_2f3.a1 + k) * (c_2f3.a2 + k) * x / ((c_2f3.b1 + k) * (c_2f3.b2 + k) * (c_2f3.b3 + k) * (k + 1));

        ++k;

        return term;
}
