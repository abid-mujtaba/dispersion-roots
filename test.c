#include <stdio.h>
#include <mpfr.h>

#include "roots.h"
#include "constants.h"
#include "dispersion.h"


#define SIZE 1


void constants();
void test1();
void test2();
void test3();
void test4();


int main(void)
{
    constants();
    test4();

    printf("\n\n");

    return 0;
}


void constants()
{
    printf("\nKAPPA_C = %.17g", KAPPA_C);
    printf("\nKAPPA_H = %.17g", KAPPA_H);
    printf("\nLAMBDA_C = %.17g", LAMBDA_C);
    printf("\nLAMBDA_H = %.17g", LAMBDA_H);
    printf("\nN0H_BY_N0E = %.17g", N0H_BY_N0E);
    printf("\nTH_BY_TC = %.17g\n", TH_BY_TC);
}


void test1()
{
    // Print samples of Disperstion function.
    for (double w0 = 1; w0 < 8; ++w0)
    {
        printf("\n");

        for (double dw = 0.1; dw < 1; dw += 0.4)
        {
            double w = w0 + dw;

            for (int k = 0; k < 40; k += 10)
            {
                    printf("\nD(%2d, %.1f) = %.17g", k, w, D(k,w));
            }
        }
    }
}


void test2()
{
    const double w = 2.25;

    printf("\nD(0, %.2f) = %.17g", w, D(0,w));
}


void test3()
{
    for (double w0 = 1; w0 < 8; ++w0, printf("\n"))
        for (double dw = 0.05; dw <= 0.95; dw += 0.05)
            printf("\nD(0, %.2f) = %.17g", w0 + dw, D(0,w0 + dw));
}


void test4()
{
    double root;
    const double k_perp = 0;
    const double omega_min = 2.05;
    const double omega_max = 2.95;

    int num = find_omega_root(k_perp, omega_min, omega_max, & root, 0.5 * (omega_min + omega_max));

    printf("\nNum of roots (between omega = %.2f and %.2f at k_perp = %.1f) = %d", omega_min, omega_max, k_perp, num);

    if (num)
        printf("\nRoot = %.17g", root);

int find_omega_root(const double k_perp, const double lo, const double hi, double * root, const double guess);
}
