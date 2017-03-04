/*
 * Implement the sine function in the range 0 to pi using the power series and
 * analytical continuation.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "sine.h"

int power_of_minus_one(int n);

double sine(double x)
{
        int n, p;
        double result = 0;

        for (n = 0; n < MAX_TERMS; ++n)
        {
                p = 2 * n + 1;

                result += power_of_minus_one(n) * gsl_pow_int(x, p) / gsl_sf_fact(p);
        }

        return result;
}

int power_of_minus_one(int n)
{
        if (n % 2 == 0)
                return 1;

        return -1;
}
