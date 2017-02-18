#include <stdio.h>
#include "functions.h"

int main(void)
{
        double x = 1.0;

        printf("\nGamma_1(%.1f) = %.18e\n", x, Gamma_n(1, x));

        return 0;
}
