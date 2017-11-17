#include <stdio.h>
#include <mpfr.h>

#include "roots.h"
#include "constants.h"
#include "dispersion.h"


#define SIZE 1


void test1();


int main(void)
{
        test1();

        printf("\n\n");

        return 0;
}


void test1()
{
        printf("\nKAPPA_C = %.17g", KAPPA_C);
        printf("\nKAPPA_H = %.17g", KAPPA_H);
        printf("\nLAMBDA_C = %.17g", LAMBDA_C);
        printf("\nLAMBDA_H = %.17g", LAMBDA_H);
        printf("\nN0H_BY_N0E = %.17g", N0H_BY_N0E);
        printf("\nTH_BY_TC = %.17g\n", TH_BY_TC);

        // Print samples of Disperstion function.
        for (double w0 = 1; w0 < 8; ++w0)
        {
                for (double dw = 0.1; dw < 1; dw += 0.4)
                {
                        double w = w0 + dw;

                        for (int k = 1; k < 40; k += 5)
                        {
                                printf("\nD(%2d, %.1f) = %.17g", k, w, D(k,w));
                        }
                }
        }
}
