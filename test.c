#include <stdio.h>
#include "functions.h"

int main(void)
{
        double xs[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
        double Gs[5];

        Gamma_n_array(1, xs, Gs, 5);

        printf("\n");

        int i;
        for (i = 0; i < 5; i++)
                printf("\t%.9e", Gs[i]);

        printf("\n\n");

        return 0;
}
