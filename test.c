#include <stdio.h>
#include "functions.h"

int main(void)
{
        int i;

        for (i = 0; i < 10; i++)
                printf("\nD(%.1f) = %.5f", i * 0.25, D(i * 0.25));

        printf("\n\n");

        return 0;
}
