#include <stdio.h>
#include "functions.h"

int main(void)
{
        double x = 5.0;
        double y = J0(x);

        printf("J0(%g) = %.18e\n", x, y);

        return 0;
}
