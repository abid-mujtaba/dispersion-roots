#include <stdio.h>
#include "functions.h"

int main(void)
{
        double x = 5.0;
        double y = J0(x);

        printf("J0(%g) = %.18e\n", x, y);

        double a[4] = {1.1, 2.2, 3.3, 4.4};

        printf("Sum = %.2f\n", sum(a, 4));

        return 0;
}
