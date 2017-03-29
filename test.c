#include <stdio.h>
#include "roots.h"


int main(void)
{
        double omega = 1.85;

        printf("\nRoot at omega = %.2f is %f", omega, find_k_perp_root(omega, 1e-3, 100));

        printf("\n\n");

        return 0;
}
