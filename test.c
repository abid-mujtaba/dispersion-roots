#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        double k_perp = 3.1;
        double omega;

        omega = 6.5 - 1e-10;
        printf("\n\nD(%.1f, %.1f) = %f", k_perp, omega, D(k_perp, omega));

        omega = 6.5 + 1e-10;
        printf("\n\nD(%.1f, %.1f) = %f", k_perp, omega, D(k_perp, omega));

        printf("\n\n");

        return 0;
}
