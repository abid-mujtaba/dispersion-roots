#include <stdio.h>
#include <pthread.h>
#include "roots.h"

// #define SIZE 100
// #define OMEGA_MAX 8
#define SIZE 3
#define OMEGA_MAX 3


// Define a struct for carrying data to and from the threads
typedef struct _thread_data {
    int start;
    double omegas[SIZE];
    double k_perps[SIZE];
    int length;
} thread_data;


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


        pthread_t threads[OMEGA_MAX - 1];       // An array containing threads
        thread_data datas[OMEGA_MAX - 1];        


        for (int i = 0; i < OMEGA_MAX - 1; ++i)
        {
            datas[i].start = i + 1;

            // pthread_create creates and executes the thread. Since 'start' is an int array 'start + i' is the pointer to its i-th element. The same is true for 'threads'
            // If the thread creation fails the function returns a non-zero value which we check for, print an error message and exit the program.
            if (pthread_create(threads + i, NULL, thread_find_omega_roots_array, datas + i))
            {
                fprintf(stderr, "Error creating thread # %d.\n", i);
                return 1;           // Failure exit code
            }
        }


        // We sequentially join with all the threads waiting for all of them to finish
        for (int i = 0; i < OMEGA_MAX - 1; ++i)
        {
            if (pthread_join(threads[i], NULL))
            {
                fprintf(stderr, "Error joining thread # %d.\n", i);
                return 2;
           }
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

    thread_data * data = (thread_data *) param;

    int N = find_omega_roots_array(data->start, k_perps, omegas, SIZE);

    for (int k = 0; k < N; ++k)
            printf("\n%d, %.1f, %.17g", data->start, k_perps[k], omegas[k]);

    pthread_exit(NULL);
}
