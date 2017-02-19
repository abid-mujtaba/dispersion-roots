#include <stdio.h>
#include "functions.h"

int main(void)
{
        double xs[5] = {0.0, 0.5, 1.5, 2.5};

        int i;

        for (i = 0; i < 4; i++)
                printf("\nSummand_1(%.1f) = %.5f", xs[i], Summand_n(1, xs[i]));

        printf("\n\n");

        return 0;
}
