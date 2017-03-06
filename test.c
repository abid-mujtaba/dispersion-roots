#include <stdio.h>
#include "hypergeom.h"


int main(void)
{
        printf("\n1F2(5; 7, 4; 1.5) = %.16f", hyp1F2(5, 7, 4, 1.5));
        printf("\n\n2F3(1, 0.5; -3.5, 2.5, -0.5; 1.5) = %.16f", hyp2F3(1, 0.5, -3.5, 2.5, -0.5, 1.5));

        printf("\n\n");

        return 0;
}
