#include <stdio.h>
#include <math.h>
#include "functions.h"
#include "hypergeom.h"


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
        //
        // printf("\n\n");
        //
        // for (i = 0; i < 20; ++i)
        // {
        //         double w = i / 20.0;
        //
        //         printf("\nD(15, %.2f) = %.5e", w, D(15, w));
        // }

        double k = 6.8, w = 0.5;
        printf("\n\nD(%.1f, %.1f) = %.20Lf", k, w, D(k, w));

        printf("\n\n");

        return 0;
}
