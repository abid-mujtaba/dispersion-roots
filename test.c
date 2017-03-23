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


        double k = 4, w = 0.5;
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

        // struct coeffs_2f3 c_2f3;
        //
        // c_2f3.a1 = 1;
        // c_2f3.a2 = 0.5;
        // c_2f3.b1 = 0.5 - KAPPA_H;
        // c_2f3.b2 = 1 + (0.5 / OMEGA_CH);
        // c_2f3.b3 = 1 - (0.5 / OMEGA_CH);
        //
        // double steps_2F3[NUM_STEPS][MAX_TERMS];
        // populate_steps_2F3(c_2f3, steps_2F3);
        //
        // int i;
        // for (i = 0; i < 10; ++i)
        //         printf("\n%d - 2F3 = %g", i, steps_2F3[1][i]);


        printf("\n\n");

        return 0;
}
