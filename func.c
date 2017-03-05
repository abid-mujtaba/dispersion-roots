/*
 * Implement the sine function in the range 0 to pi using the power series and
 * analytical continuation.
 */

#include <stdio.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "func.h"

int power_of_minus_one(int n);
double a(int k);                // Calculate the coefficient a_k for (x - 3/2)^k expansion
double ln_3_2(double x);

double ln(double x)
{
        if (x > 2)
                return ln_3_2(x);

        int n;
        double result = 0;
        double y = x - 1;               // Because of the way power series of ln is defined

        for (n = 1; n < MAX_TERMS; ++n)
        {
                result += power_of_minus_one(n + 1) * gsl_pow_int(y, n) / n;
        }

        return result;
}


double ln_3_2(double x)
{
        int k;
        double result = 0;
        double y = x - 1.5;

        printf("\n\n");

        for (k = 0; k < MAX_TERMS; ++k)
                result += a(k) * gsl_pow_int(y, k);

        printf("\n");

        return result;
}


int power_of_minus_one(int n)
{
        if (n % 2 == 0)
                return 1;

        return -1;
}


double a(int k)
{
        setbuf(stdout, NULL);           // Disable buffering on the stdout stream

        double result = 0;
        double term;
        int n;
        int start = k ? k : 1;          // If k = 0 then start should be 1 for calculating a_0 since the zeroeth term is not part of
                                        // the expansion of a_0 and it blows up (-inf)

        for (n = start; n < start + MAX_TERMS; ++n)
        {
                // result += power_of_minus_one(n + 1) * gsl_sf_choose(n, k) / (n * gsl_pow_int(2, n - k));
                term = power_of_minus_one(n + 1) * gsl_sf_choose(n, k) / (n * gsl_pow_int(2, n - k));

                if (k == 25)
                {
                        // We use \r to print on the same line repeatedly
                        printf("\rb(%d) = %-20.9f", n, term);
                        // fflush(stdout);         // Use fflush to force the buffer to flush. Not needed if using setbuf(stdout, NULL);
                        usleep(5e4);
                }

                result += term;
        }

        return result;
}
