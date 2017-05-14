#include <stdio.h>
#include "roots.h"
#include "constants.h"
#include <mpfr.h>


void test1();


int main(void)
{
        // test1();
        mpfr_t x;
        mpfr_init(x);

        mpfr_set_d(x, 7.42, RND);
        mpfr_gamma(x, x, RND);
        mpfr_printf("\nx = %.17RG", x);

        char * res;
        mpfr_exp_t * expptr;

        res = mpfr_get_str((char *) 0, expptr, 10, 0, x, RND);
        printf("\nstring = %s", res);
        mpfr_free_str(res);

        mpfr_clear(x);


        printf("\n\n");

        return 0;
}


void test1()
{
        const int size = 20;
        const int initial = 6;
        double k_perps[size];
        double omegas[size];
        int num;

        num = find_omega_roots_array(initial, k_perps, omegas, size);

        for (int i = 0; i < num; ++i)
                printf("\nD root at k_perp = %.2f  -> %.17g", k_perps[i], omegas[i]);
}
