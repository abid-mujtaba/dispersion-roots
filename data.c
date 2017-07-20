#include <stdio.h>
#include <pthread.h>
#include "roots.h"

// #define SIZE 100
// #define OMEGA_MAX 8
#define SIZE 3
#define OMEGA_MAX 3


void * thread_find_omega_roots_array(void * param);


int main(void)
{
        FILE * fout = fopen("data.csv", "w");
        fprintf(fout, "seq,k_perp,omega");

        /*double omegas[SIZE];*/
        /*double k_perps[SIZE];*/

        /*for (int i = 1; i < OMEGA_MAX; ++i)*/
        /*{*/
                /*int N = find_omega_roots_array(i, k_perps, omegas, SIZE);*/

                /*for (int k = 0; k < N; ++k)*/
                        /*fprintf(fout, "\n%d,%.1f,%.17g", i, k_perps[k], omegas[k]);*/
        /*}*/


        pthread_t threads[OMEGA_MAX - 1];
        int start[OMEGA_MAX - 1];

        for (int i = 0; i < OMEGA_MAX - 1; ++i)
        {
            start[i] = i + 1;
            pthread_create(& threads[i], NULL, thread_find_omega_roots_array, start + i);
        }


        for (int i = 0; i < OMEGA_MAX - 1; ++i)
        {
            pthread_join(threads[i], NULL);
        }



        fprintf(fout, "\n");
        fclose(fout);

        printf("\n\n");

        return 0;
}


void * thread_find_omega_roots_array(void * param)
{
    double omegas[SIZE];
    double k_perps[SIZE];

    int * omega_start_ptr = (int *) param;

    int N = find_omega_roots_array(* omega_start_ptr, k_perps, omegas, SIZE);

    for (int k = 0; k < N; ++k)
            printf("\n%d, %.1f, %.17g", * omega_start_ptr, k_perps[k], omegas[k]);

    return NULL;
}
