/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <math.h>
#include "hypergeom.h"


double hyp1F2(const double a1, const double b1, const double b2, const double x)
{
        int k;


        double result = 0;
        double term = 1;                // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE

        // When (the absolute value of) 'term' becomes less than TOLERANCE the sum stops changing rapidly and we truncate it
        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                result += term;

                // Based on the definition of 1F2 each term differs from the previous by multiplicative factors that have to do with the Pochhammer symbols, power of x and the factorial in the denominator
                term *= (a1 + k) * x / ((b1 + k) * (b2 + k) * (k + 1));
        }

        return result;
}


double hyp2F3(const double a1, const double a2, const double b1, const double b2, const double b3, const double x)
{
        int k;

        double result = 0;
        double term = 1;

        for (k = 0; (k < MAX_TERMS) & (fabs(term) > TOLERANCE); ++k)
        {
                result += term;

                term *= (a1 + k) * (a2 + k) * x / ((b1 + k) * (b2 + k) * (b3 + k) * (k + 1));
        }

        return result;
}
