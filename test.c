#include <stdio.h>
#include "func.h"


int main(void)
{
        double x1 = 1.25;
        double x2 = 1.5;
        double x3 = 2.25;

        printf("\nln(%.2f) = %.9f", x1, ln(x1));
        printf("\nln(%.2f) = %.9f", x2, ln(x2));
        printf("\nln(%.2f) = %.9f", x3, ln(x3));
        printf("\n\n");

        return 0;
}
