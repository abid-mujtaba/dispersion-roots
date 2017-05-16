#include <stdio.h>
#include "roots.h"
#include "constants.h"
#include <mpfr.h>


void test1();
void test2();


int main(void)
{
        // test1();

        struct Constants cc;
        get_constants_c(& cc);
        printf("\nRHO_C = %f", cc.rho);

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


// gamma(7.42) in 1024 bit precision
#define STR "1.602117949820571510274107883952300756058764977396274970355879617540504416806075178095409660181443828643423198484563617851743062024608885523111927964457849198546580997114798115486670867610256520981302272757307249459052852987114940791087119415189679469472082755100956878053655788810840623971816764448986572307187e3"


void test2()
{
        mpfr_set_default_prec(512);

        mpfr_t x, y;
        mpfr_init(x);
        mpfr_init(y);

        mpfr_set_d(x, 7.42, RND);
        mpfr_gamma(x, x, RND);
        mpfr_printf("\nx = %.17RG\n", x);

        // Print the value of gamma(7.42) to default precision
        printf("\nGamma(7.42) as string = ");
        mpfr_out_str((FILE *) 0, 10, 0, x, RND);                // Null FILE pointer means it will be output to \n

        mpfr_set_str(y, STR, 10, RND);          // Set value of y using the string defined above
        mpfr_printf("\n\ny = %.17RG\n", y);

        mpfr_sub(x, x, y, RND);
        mpfr_printf("\nx - y = %.17RG", x);     // Subtract x and y to show that they are equal

        mpfr_clear(x);
        mpfr_clear(y);
}
