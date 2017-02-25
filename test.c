#include <stdio.h>
#include <math.h>
#include "functions.h"

int main(void)
{
        double a = 1.7;
        double b = -1.3;

        int x = signbit(a);
        int y = signbit(b);

        printf("%.2f XOR %.2f = %d XOR %d = %d", a, b, x, y, x ^ y);

        if (x ^ y)
                printf("\nTrue");
        else
                printf("\nFalse");

        printf("\n\n");

        return 0;
}
