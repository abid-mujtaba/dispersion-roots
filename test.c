#include <stdio.h>
#include "functions.h"

int main(void)
{
        double x = 5.0;
        double y = J0(x);

        printf("\nJ0(%g) = %.18e\n", x, y);

        double a[4] = {1.1, 2.2, 3.3, 4.4};
        printf("Sum = %.2f\n", sum(a, 4));

        double b[4];
        triple(a, b, 4);
        printf("\nTripled:");

        int i = 0;

        for (i = 0; i < 4; i++)
                printf("\t%.2f", b[i]);

        printf("\n\n");

        double xs[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        double js[10];
        J0_array(xs, js, 10);

        for (i = 0; i < 10; i++)
                printf("  %.6f", js[i]);

        printf("\n\n");

        return 0;
}
