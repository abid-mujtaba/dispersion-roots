#include <stdio.h>
#include "functions.h"

int main(void)
{
        double xs[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        double js[10];
        J0_array(xs, js, 10);

        int i;

        for (i = 0; i < 10; i++)
                printf("  %.6f", js[i]);

        printf("\n\n");

        return 0;
}
