#include <stdio.h>
#include "roots.h"

#define SIZE 200


int main(void)
{
        FILE * fout = fopen("data.csv", "w");
        fprintf(fout, "seq,k_perp,omega");

        double omegas[SIZE];
        double k_perps[SIZE];

        int N = find_omega_roots_array(1, k_perps, omegas, SIZE);

        for (int k = 0; k < N; ++k)
                fprintf(fout, "\n1,%.1f,%.17g", k_perps[k], omegas[k]);

        fprintf(fout, "\n");
        fclose(fout);

        return 0;
}
