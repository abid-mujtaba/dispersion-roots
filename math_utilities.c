// Implement functions for common math utilities.

#include <mpfr.h>
#include "constants.h"


int factorial(int n)
{
    if (n == 0)
        return 1;

    return n * factorial(n - 1);
}


// Pochammer symbol with negative index defined as
// (a)_-n = (a) (a - 1) ... (a - n + 1)
void neg_pochammer(mpfr_t res, int n, const mpfr_t a, mpfr_t x)
{
    mpfr_set_ui(res, 1, RND);       // (a)_0

    for (int i = 1; i <= n; ++i)
    {
        mpfr_sub_ui(x, a, -(i - 1), RND);
        mpfr_mul(res, res, x, RND);         // res *= (a - (i - 1))
    }
}
