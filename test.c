#include <stdio.h>
#include <mpfr.h>

#include "roots.h"
#include "constants.h"
#include "dispersion.h"


#define SIZE 1


void test1();
void test2();
void test3();


int main(void)
{
        test3();

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
    const int initial = 5;

    double samples[SIZE] = {0.0};

    double k_perps[SIZE];
    double omegas[SIZE];

    int num = find_omega_roots_array(initial, samples, k_perps, omegas, SIZE);

    for (int i = 0; i < num; ++i)
        printf("\nRoot at k_perp = %.2f  ->  %.2f", k_perps[i], omegas[i]);
}


void test3()
{
    const double k_perp = 0;
    double omega;

    for (int i = 1; i < 8; ++i)
    {
        printf("\n");

        for (int j = 0; j < 10; ++j)
        {
            omega = i + j / 10.0;
            printf("\nD(%.2f, %.2f) = %.5f", k_perp, omega, D(k_perp, omega));
        }
    }
}
