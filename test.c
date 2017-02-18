#include <stdio.h>
#include "functions.h"

int main(void)
{
        int n = 0;

        double xs[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        double Ins[10];
        In_array(n, xs, Ins, 10);

        int i;

        for (i = 0; i < 10; i++)
                printf("  %.6f", Ins[i]);

        printf("\n\n");

        return 0;
}
