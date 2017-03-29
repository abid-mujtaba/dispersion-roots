#include <stdio.h>
#include "functions.h"


int main(void)
{
        double k = 0;

        while (k <= 100)
        {
            double k_perp = k;
            double omega = 1.85;

            k += 10;

            printf("\n\nk_perp = %f", k_perp);
            printf("\nD(%.1f, %.1f) = %f", k_perp, omega, D(k_perp, omega));
        }

        printf("\n\n");

        return 0;
}
