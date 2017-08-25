#include <stdio.h>
#include <mpfr.h>

#include "roots.h"
#include "constants.h"
#include "dispersion.h"

#include "alpha.h"


#define SIZE 1


void test1();
void test2();
void test3();
void test4();
void test5();
void test6();


int main(void)
{
        test6();

        printf("\n\n");

        return 0;
}


void test1()
{
        /*const int size = 20;*/
        /*const int initial = 6;*/
        /*double k_perps[size];*/
        /*double omegas[size];*/
        /*int num = 0;*/

        /*num = find_omega_roots_array(initial, k_perps, omegas, size);*/

        /*for (int i = 0; i < num; ++i)*/
                /*printf("\nD root at k_perp = %.2f  -> %.17g", k_perps[i], omegas[i]);*/
}


void test2()
{
    const double omega = 1.5;
    double k_perp;

    for (int i = 0; i < 10; i++)
    {
        k_perp = i * 10 + 1.0;
        printf("\nD(%.2f, %.2f) = %.9f", k_perp, omega, D(k_perp, omega));
    }
}


void test4()
{
    const double k_perp = 7.5;
    double omega;

    for (int i = 0; i < 7; ++i)
    {
        omega = 1.5 + i;
        printf("\nD(%.2f, %.3f) = %.9f", k_perp, omega, D(k_perp, omega));
    }
}



void test3()
{
    const int initial = 1;

    double samples[SIZE] = {2.6};

    double k_perps[SIZE];
    double omegas[SIZE];

    int num = find_omega_roots_array(initial, samples, k_perps, omegas, SIZE);

    for (int i = 0; i < num; ++i)
        printf("\nRoot at k_perp = %.2f  ->  %.2f", k_perps[i], omegas[i]);
}


void test5()
{
    mpfr_t kappa, result, x, y;
    mpfr_inits(kappa, result, x, y, (mpfr_ptr) 0);

    mpfr_set_d(kappa, KAPPA_H, RND);
    alpha(result, 1, LAMBDA, kappa, x, y);

    mpfr_printf("\nalpha(%d) = %RG", 1, result);

    mpfr_clears(kappa, result, x, y, (mpfr_ptr) 0);
}


void test6()
{
    double k_perp = 1;
    double omega = 1.5;

    printf("\nD(%.2f, %.2f) = %.9f", k_perp, omega, D(k_perp, omega));
}
