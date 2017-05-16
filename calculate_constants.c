/*
 * Calculate the constants derived from the ones defined in constants.h and store
 * them in derived.h.
 *
 * This process will be carried out before compilation so that these values are
 * pre-compiled in to the executable and the calculations won't have to be
 * carried out over and over again.
 */

#include <stdio.h>
#include <math.h>
#include "constants.h"


int main(void)
{
        FILE * fout = fopen("derived.h", "w");

        fprintf(fout, "#define RHO_C %.17g", RHO_H / sqrt(TH_BY_TC));            // .17g guarantees that the full double is printed
        fprintf(fout, "\n#define N0C_BY_N0E %.17g", 1 - N0H_BY_N0E);

        fprintf(fout, "\n");
        fclose(fout);
}
