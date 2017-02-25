#include <stdio.h>
#include "functions.h"

int main(void)
{
        double slices[8] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2};
        double omega[16];
        double k_perp[16];

        int i;
        int num = find_k_perp_roots_array(slices, omega, k_perp, 8);

        for (i = 0; i < num; i++)
                printf("\nRoot at %.2f = %.5f", omega[i], k_perp[i]);

        printf("\n\n");

        return 0;
}
