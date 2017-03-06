#include <stdio.h>
#include "hypergeom.h"


int main(void)
{
        printf("\n1F2(3; 5, 2; 1.5) = %.15f", hyp1F2(3, 5, 2, 1.5));
        printf("\n\n2F3(1, 0.5; -1.5, 2.5, -0.5; 1.5) = %.15f", hyp2F3(1, 0.5, -1.5, 2.5, -0.5, 1.5));

        printf("\n\n");

        return 0;
}
