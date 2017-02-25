#include <stdio.h>
#include "functions.h"

#define UPPER 4.5
#define INTERVAL 100.0
#define NUM (int) (UPPER * INTERVAL)

int main(void)
{
        double slices[NUM];
        double omega[NUM * 2], k_perp[NUM * 2];

        int i;
        for (i = 0; i < NUM; ++i)
                slices[i] = i / INTERVAL;

        int num = find_k_perp_roots_array(slices, omega, k_perp, NUM);

        printf("Total number of roots found = %d", num);

        printf("\n\n");

        return 0;
}
