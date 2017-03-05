/*
 * Implement the sine function in the range 0 to pi using the power series and
 * analytical continuation.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "func.h"

int power_of_minus_one(int n);

double ln(double x)
{
        int n;
        double result = 0;
        double y = x - 1;               // Because of the way power series of ln is defined

        for (n = 1; n < MAX_TERMS; ++n)
        {
                result += power_of_minus_one(n + 1) * gsl_pow_int(y, n) / n;
        }

        return result;
}

int power_of_minus_one(int n)
{
        if (n % 2 == 0)
                return 1;

        return -1;
}
