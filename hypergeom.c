/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include "hypergeom.h"


double hyp1F2(const double a1, const double b1, const double b2, const double x)
{
        int k;

        // Define variables that will contain the running value of the Pochhammer symbols corresponding to the a_i and b_i
        double pa1 = 1;         // (a1)_0 = 1 by the definition of the Pochhammer symbol
        double pb1 = 1;
        double pb2 = 1;

        // Define running variables for power of x and factorial
        double fact = 1;
        double pow_x = 1;

        double result = 0;

        for (k = 0; k < MAX_TERMS; ++k)
        {
                result += pa1 * pow_x / (pb1 * pb2 * fact);

                // Update all running variables
                pa1 *= a1 + k;
                pb1 *= b1 + k;
                pb2 *= b2 + k;

                fact *= (k + 1);
                pow_x *= x;
        }

        return result;
}


double hyp2F3(const double a1, const double a2, const double b1, const double b2, const double b3, const double x)
{
        int k;

        // Define variables that will contain the running value of the Pochhammer symbols corresponding to the a_i and b_i
        double pa1 = 1;         // (a1)_0 = 1 by the definition of the Pochhammer symbol
        double pa2 = 1;
        double pb1 = 1;
        double pb2 = 1;
        double pb3 = 1;

        // Define running variables for power of x and factorial
        double fact = 1;
        double pow_x = 1;

        double result = 0;

        for (k = 0; k < MAX_TERMS; ++k)
        {
                result += pa1 * pa2 * pow_x / (pb1 * pb2 * pb3 * fact);

                // Update all running variables
                pa1 *= a1 + k;
                pa2 *= a2 + k;
                pb1 *= b1 + k;
                pb2 *= b2 + k;
                pb3 *= b3 + k;

                fact *= (k + 1);
                pow_x *= x;
        }

        return result;
}
