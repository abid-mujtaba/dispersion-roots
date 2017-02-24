#include <stdio.h>
#include "functions.h"

int main(void)
{
        int i;
        struct D_params d_params;

        for (i = 0; i < 10; i++)
        {
                d_params.omega = i * 0.25;
                printf("\nD(1, %.1f) = %.5f - D_root(1, %.1f) = %.5f", i * 0.25, D(K_PERP, i * 0.25), d_params.omega, D_root(K_PERP, &d_params));
        }

        printf("\n\n");

        return 0;
}
