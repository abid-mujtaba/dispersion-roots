#include <stdio.h>
#include <math.h>
#include "functions.h"


int main(void)
{
        int i;

        for (i = 0; i < 20; ++i)
        {
                double k = i / 5.0;

                printf("\nD(%.2f, 1.5) = %.5e", k, D(k, 1.5));
        }

        printf("\n\n");

        for (i = 0; i < 20; ++i)
        {
                double w = i / 5.0;

                printf("\nD(1.5, %.2f) = %.5e", w, D(1.5, w));
        }

        printf("\n\n");

        return 0;
}
