#include <stdio.h>
#include "roots.h"
#include "functions.h"


int main(void)
{
        double slices[10];
        double omegas[10];
        double roots[10];
        int NUM = 7;

        for (int j = 0; j < NUM; ++j)
                slices[j] = j + 1.5 + 1e-10;

        int N = find_k_perp_roots_array(slices, omegas, roots, NUM);

        for (int i = 0; i < N; ++i)
                printf("\nRoot at omega = %.2f is %f", omegas[i], roots[i]);

        printf("\n\n");

        return 0;
}
