#include <stdio.h>
#include "roots.h"


int main(void)
{
        const int size = 10;
        const int initial = 6;
        double k_perps[size];
        double omegas[size];
        int num;

        num = find_omega_roots_array(initial, k_perps, omegas, size);

        for (int i = 0; i < num; ++i)
                printf("\nD root at k_perp = %.2f  -> %.17g", k_perps[i], omegas[i]);

        printf("\n\n");

        return 0;
}
