#include <stdio.h>
#include "roots.h"


int main(void)
{
        double lo = 2 + 1e-6;
        double hi = 3 - 1e-6;
        double k_perp = 1.6;
        double root;

        find_omega_root(k_perp, lo, hi, & root);
        printf("\nRoot of D between %.2f and %.2f at k_perp = %.2f  ->  %.17g", lo, hi, k_perp, root);

        printf("\n\n");

        return 0;
}
