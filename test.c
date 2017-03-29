#include <stdio.h>
#include "roots.h"


int main(void)
{
        double omega = 1.05;

        while (omega < 2)
        {
                printf("\nRoot at omega = %.2f is %f", omega, find_k_perp_root(omega, 1e-3, 100));

                omega += 0.05;
        }

        printf("\n\n");

        return 0;
}
