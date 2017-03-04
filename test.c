#include <stdio.h>
#include "roots.h"

#include <gsl/gsl_sf_gamma.h>

#define UPPER 4.5
#define INTERVAL 100.0
#define NUM (int) (UPPER * INTERVAL)

int main(void)
{
        // double slices[NUM];
        // double omega[NUM * 2], k_perp[NUM * 2];
        //
        // int i;
        // for (i = 0; i < NUM; ++i)
        //         slices[i] = i / INTERVAL;
        //
        // int num = find_k_perp_roots_array(slices, omega, k_perp, NUM);
        //
        // printf("Total number of roots found = %d", num);


        int i;
        const int s = 5;
        double xs[5] = {-1.5, -0.5, 0.5, 1, 2};
        double ys[s];

        // for (i = 0; i < s; ++i)
        // {
        //         xs[i] = -(i+1) / 5.0;
        //         printf("\nx = %.2f", xs[i]);
        //         printf("\nGamma(%.2f) = %.2f\n", xs[i], gsl_sf_gamma(xs[i]));
        // }

        Gamma_array(xs, ys, s);

        for (i = 0; i < s; ++i)
                printf("\nGamma(%.2f) = %.2f", xs[i], ys[i]);

        printf("\n\n");

        return 0;
}
