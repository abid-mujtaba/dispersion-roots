#include <stdio.h>
#include <math.h>
#include "functions.h"
#include "hypergeom.h"
#include "constants.h"


int main(void)
{
        // int i;
        //
        // for (i = 0; i < 100; ++i)
        // {
        //         double k = i / 5.0;
        //
        //         printf("\nD(%.2f, 0.5) = %.5e", k, D(k, 0.5));
        // }

        // for (i = 0; i < 20; ++i)
        // {
        //         double w = i / 20.0;
        //
        //         printf("\nD(15, %.2f) = %.5e", w, D(15, w));
        // }

        double k = 6.8, w = 0.5;
        printf("\n\nD(%.1f, %.1f) = %.20f", k, w, D(k, w));


        // struct coeffs_1f2 c_1f2;
        //
        // c_1f2.a1 = KAPPA_H + 1;
        // c_1f2.b1 = KAPPA_H + 1.5 + (0.5 / OMEGA_CH);
        // c_1f2.b2 = KAPPA_H + 1.5 - (0.5 / OMEGA_CH);
        //
        // // We create and populate the intermediate step values in the hypergeometric functions
        // double steps_1F2[NUM_STEPS][MAX_TERMS];
        // populate_steps_1F2(c_1f2, steps_1F2);

        printf("\n\n");

        return 0;
}
