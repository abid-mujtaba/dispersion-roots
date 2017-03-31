#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        int NUM = 20;

        double slices[NUM];
        double omegas[NUM];
        double roots[NUM];

        for (int j = 0; j < NUM; ++j)
                slices[j] = j * 0.5 + 1.1 + 1e-10;

        int N = find_k_perp_roots_array(slices, omegas, roots, NUM);

        for (int i = 0; i < N; ++i)
                printf("\n%2d. - Root at omega = %.2f is %f", i, omegas[i], roots[i]);

        printf("\n\n");

        return 0;
}
