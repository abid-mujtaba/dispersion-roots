#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        const int size = 200;
        double omegas[size];
        double k_perps[size];

        int N = find_omega_roots_array(1, k_perps, omegas, size);

        for (int k = 0; k < N; ++k)
                printf("\nRoot at k_perp = %.2f  =  %f", k_perps[k], omegas[k]);

        printf("\n\n");

        return 0;
}
