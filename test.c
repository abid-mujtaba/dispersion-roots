#include <stdio.h>
#include "constants.h"
#include "functions.h"


int main(void)
{
        printf("\nRHO_C = %.20f", RHO_C);
        printf("\nRHO_H = %d", RHO_H);
        printf("\nLAMBDA_KAPPA_C_p2 = %.20f", LAMBDA_KAPPA_C_p2);
        printf("\nLAMBDA_KAPPA_H_p2 = %.20f", LAMBDA_KAPPA_H_p2);
        printf("\n\nspecie_c(1.5, 1.5) = %.20f", specie_c(1, 1.5));
        printf("\n\nspecie_h(1.5, 1.5) = %.20f", specie_h(1, 1.5));

        printf("\n\n");

        return 0;
}
