#include <stdio.h>
#include "functions.h"


// Declare the function D to refer to the Henning Dispersion relation/function
double (* D)(double, double) = D_Henning;


int main(void)
{
        double k_perp = 3.1;
        double om1 = 6.2;
        double om2 = 6.5;

        printf("\nD(%.1f, %.1f) = %.17g", k_perp, om1, D(k_perp, om1));
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, om2, D(k_perp, om2));

        printf("\n\n");

        return 0;
}
