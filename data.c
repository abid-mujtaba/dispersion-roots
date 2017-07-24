#include <stdio.h>
#include <pthread.h>
#include "roots.h"

#define SIZE 100
#define OMEGA_MAX 8

// Since we are studying intervals of omega starting at 1 the number of threads is one less than OMEGA_MAX
#define NUM_THREADS (OMEGA_MAX - 1)


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


        pthread_t threads[NUM_THREADS];       // An array containing threads
        thread_data datas[NUM_THREADS];        


        for (int i = 0; i < NUM_THREADS; ++i)
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
        for (int i = 0; i < NUM_THREADS; ++i)
        {
            if (pthread_join(threads[i], NULL))
            {
                fprintf(stderr, "Error joining thread # %d.\n", i);
                return 2;
            }
        }


        // Store calculated data (which is not stored in the 'datas' array) to the file
        for (int i = 0; i < NUM_THREADS; ++i)
        {
            thread_data data = datas[i];

            for (int k = 0; k < data.length; ++k)
                fprintf(fout, "\n%d,%.1f,%.17g", data.start, data.k_perps[k], data.omegas[k]);
        }



        fprintf(fout, "\n");
        fclose(fout);

        printf("\n");

        return 0;
}


void * thread_find_omega_roots_array(void * param)
{
    thread_data * data = (thread_data *) param;

    int N = find_omega_roots_array(data->start, data->k_perps, data->omegas, SIZE);
    data->length = N;        // Store the number of calculated roots in the struct

    return NULL;
}
