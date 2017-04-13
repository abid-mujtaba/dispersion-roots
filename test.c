#include <stdio.h>
#include "dispersion.h"
#include "henning.h"


int main(void)
{
        double k_perp = 3.1;
        double omega = 6.2;

        printf("\nD_Henning(%.1f, %.1f) = %.17g", k_perp, omega, D_Henning(k_perp, omega));
        printf("\nD(%.1f, %.1f) = %.17g", k_perp, omega, D(k_perp, omega));

        printf("\n\n");

        return 0;
}
