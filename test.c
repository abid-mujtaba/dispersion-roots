#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        double k_perp;

        for (k_perp = 0.25; k_perp < 4; k_perp += 0.25)
        {
                double omega = find_omega_root(k_perp, 1 + 1e-10, 2 - 1e-10);
                printf("\nOmega root at k_perp = %.2f  =  %f", k_perp, omega);
        }

        printf("\n\n");

        return 0;
}
