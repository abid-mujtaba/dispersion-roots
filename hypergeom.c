/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <math.h>
#include "hypergeom.h"


int terms_2F3(const struct coeffs_2f3, const double x, long double terms[]);
int terms_1F2(const double coeff, const struct coeffs_1f2, const double x, long double terms[], const int start);


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

        long double terms[2 * MAX_TERMS];

        int size = terms_2F3(c_2f3, x, terms);
        int mid = size;
        size = terms_1F2(coeff, c_1f2, x, terms, size);

        // long double result = coeff * hyp1F2(c_1f2.a1, c_1f2.b1, c_1f2.b2, x);
        long double result = 0;

        for (k = 0; k < size; ++k)
                result += terms[k];

        long double r_2f3 = 0;
        for (k = 0; k < mid; ++k)
                r_2f3 += terms[k];

        printf("\n2F3 = %.8Le", hyp2F3(c_2f3.a1, c_2f3.a2, c_2f3.b1, c_2f3.b2, c_2f3.b3, x));
        printf("\narray = %.8Le", r_2f3);
        printf("\ndiff = %.8Le\n", r_2f3 + hyp2F3(c_2f3.a1, c_2f3.a2, c_2f3.b1, c_2f3.b2, c_2f3.b3, x));


        long double r_1f2 = 0;
        for (k = mid; k < size; ++k)
                r_1f2 += terms[k];

        printf("\ncoeff * 1F2 = %.8Le", coeff * hyp1F2(c_1f2.a1, c_1f2.b1, c_1f2.b2, x));
        printf("\narray = %.8Le", r_1f2);
        printf("\ndiff = %.8Le", r_1f2 - coeff * hyp1F2(c_1f2.a1, c_1f2.b1, c_1f2.b2, x));

        return result;
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

/*
 * The same array is passed which has already been populated by terms_2F3.
 * The int 'start' tells you where to start adding terms after the earlier population.
 */
int terms_1F2(const double coeff, const struct coeffs_1f2 c, const double x, long double terms[], const int start)
{
        int k;

        long double term = coeff;

        for (k = 0; ((k < MAX_TERMS) & (fabs(term) > TOLERANCE)); ++k)
        {
                terms[start + k] = term;

                term *= (c.a1 + k) * x / ((c.b1 + k) * (c.b2 + k) * (k + 1));
        }

        return k + start;
}
