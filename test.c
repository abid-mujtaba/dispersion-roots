#include <stdio.h>
#include "functions.h"

int main(void)
{
        // double slices[8] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
        // double omega[8];
        // double k_perp[8];
        //
        // int i;
        // int num = find_k_perp_roots(slices, omega, k_perp, 8);
        //
        // for (i = 0; i < num; i++)
        //         printf("\nRoot at %.2f = %.5f", omega[i], k_perp[i]);

        find_k_perp_root(3.1);

        printf("\n\n");

        return 0;
}
