#include <stdio.h>
#include "dispersion.h"


// Declare the function D to refer to the Henning Dispersion relation/function
double (* Df)(double, double) = D;


int main(void)
{
        Df(3.1, 6.2);                // Call the Df function. All tests are defined inside it while we develop

        printf("\n\n");

        return 0;
}
