#include <stdio.h>
#include "functions.h"


int main(void)
{
        int i;

        for (i = 0; i < 20; ++i)
                printf("\nspecie_h(%.2f, 1.5) = %.20f", i / 5.0, specie_h(i / 5.0, i / 5.0));

        printf("\n\n");

        return 0;
}
