#include <stdio.h>
#include "roots.h"
#include "constants.h"
#include <mpfr.h>

#define SIZE 1


void test1();
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
