#include <stdio.h>
#include "functions.h"

int main(void)
{
        printf("\nGamma_1(0) = %.5f", Gamma_n(1, 0));
        printf("\nGamma_1(0) / 0^2 = %.5f", Gamma_n(1, 0) / (0 * 0));
        printf("\nGamma_1_by_x2(0) = %.5f", Gamma_n_by_x2(1, 0));
        printf("\nSummand_1(0, 1.5) = %.5f", Summand_n(1, 0, 1.5));
        printf("\nD(0, 1.5) = %.5f", D(0, 1.5));

        printf("\n\n");

        return 0;
}
