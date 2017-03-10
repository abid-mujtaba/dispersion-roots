#include <stdio.h>
#include <math.h>
#include "functions.h"


int main(void)
{
        int i;

        // for (i = 0; i < 20; ++i)
        //         printf("\nspecie_h(%.2f, 1.5) = %.20f", i / 5.0, specie_h(i / 5.0, 1.5));

        for (i = 0; i < 20; ++i)
        {
                double k = pow(10, -i);
                printf("\nspecie(%.1e, 1.5) = %.10e", k, specie_h(k, 1.5));
        }

        // printf("\nspecie_h(0, 1.5) = %.20f", specie_h(1e-9, 1.5));
        printf("\n\n");

        return 0;
}
