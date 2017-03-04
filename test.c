#include <stdio.h>
#include "sine.h"


int main(void)
{
        double x1 = 0.5;
        double x2 = 3;

        printf("\nSin(%.2f) = %.9f", x1, sine(x1));
        printf("\nSin(%.2f) = %.9f", x2, sine(x2));
        printf("\n\n");

        return 0;
}
