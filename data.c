#include <stdio.h>
#include "roots.h"

#define SIZE 200


int main(void)
{
        FILE * fout = fopen("data.csv", "w");
        fprintf(fout, "seq,k_perp,omega");

        double omegas[SIZE];
        double k_perps[SIZE];

        for (int i = 1; i < 8; ++i)
        {
                int N = find_omega_roots_array(i, k_perps, omegas, SIZE);

                for (int k = 0; k < N; ++k)
                        fprintf(fout, "\n%d,%.1f,%.17g", i, k_perps[k], omegas[k]);
        }

        fprintf(fout, "\n");
        fclose(fout);

        return 0;
}
