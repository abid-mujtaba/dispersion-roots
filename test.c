#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        // double slices[10];
        // double omegas[10];
        // double roots[10];
        // int NUM = 7;
        //
        // for (int j = 0; j < NUM; ++j)
        //         slices[j] = j + 1.5;
        //
        // int N = find_k_perp_roots_array(slices, omegas, roots, NUM);
        //
        // for (int i = 0; i < N; ++i)
        //         printf("\nRoot at omega = %.2f is %f", omegas[i], roots[i]);

        // double k = 0;
        // double omega = 5 + 1e-20;
        //
        // printf("\nD(%.2f, %.2f) = %f\n", k, omega, D(k, omega));
        //
        // k = 1e-20;
        // printf("\nD(%.2f, %.2f) = %f", k, omega, D(k, omega));

        double k = 0;

        double omega = 5 + 1e-10;

        for (int i = 0; i < 10; ++i)
                printf("\nD(%.2f, %.2f) = %f", k, omega + i * 0.1, D(k, omega + i * 0.1));

        printf("\n\n");

        return 0;
}
