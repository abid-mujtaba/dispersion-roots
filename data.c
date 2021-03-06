#include <stdio.h>
#include <pthread.h>
#include "roots.h"
#include "constants.h"

// Note: Calculate the following carefully based on the delta and split of k_perp values as we move from 0 to K_PERP_MAX
#define MAX_K_PERP_SAMPLES 250             // Number of samples of k_perp at which the root will be calculated.

// Starting value of omega by omega_ce
#define OMEGA_MIN 1
#define OMEGA_MAX OMEGA_MIN + 7

// The number of threads required is equal to the number of bands (omega integer gaps) we are calculating over
#define NUM_THREADS (OMEGA_MAX - OMEGA_MIN + 1)

// Declare the file names for storing results
#define VALUEFILE "data/value-" PLOT ".json"
#define DATAFILE "data/data-" PLOT ".csv"


// Define a struct for carrying data to and from the threads
typedef struct _thread_data {
    int start;
    double omegas[MAX_K_PERP_SAMPLES];
    double k_perps[MAX_K_PERP_SAMPLES];
    int length;
} thread_data;


void write_values();
void * thread_find_omega_roots_array(void * param);


int main(void)
{
        // Write various constant values to file
        write_values();

        FILE * fout = fopen(DATAFILE, "w");
        fprintf(fout, "seq,k_perp,omega");

        pthread_t threads[NUM_THREADS];       // An array containing threads
        thread_data datas[NUM_THREADS];


        for (int i = 0; i < NUM_THREADS; ++i)
        {
            datas[i].start = i + OMEGA_MIN;

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

            if (data.length > MAX_K_PERP_SAMPLES)
                fprintf(stderr, "\nWarning - Count = %d > MAX_K_PERP_SAMPLES = %d for initial omega = %d", data.length, MAX_K_PERP_SAMPLES, OMEGA_MIN + i);


            for (int k = 0; k < data.length; ++k)
                fprintf(fout, "\n%d,%.17f,%.17g", data.start, data.k_perps[k], data.omegas[k]);
        }



        fprintf(fout, "\n");
        fclose(fout);

        printf("\n");

        return 0;
}


void * thread_find_omega_roots_array(void * param)
{
    thread_data * data = (thread_data *) param;

    // We must construct the k_perp sample values which we split in to two halves to get better precision for smaller values and not waste processing on the mundane smoother slope larger values

    double k_perp_samples[MAX_K_PERP_SAMPLES];
    double delta = 0.005 * K_PERP_MAX;             // Upto 40% of max
    double sample = 0;
    int count = 0;

    while (sample < 40 && sample < K_PERP_MAX)
    {
        k_perp_samples[count++] = sample;
        sample += delta;
    }

    delta = 0.02 * K_PERP_MAX;      // After 40

    while (sample < K_PERP_MAX)
    {
        k_perp_samples[count++] = sample;
        sample += delta;
    }

    int N = find_omega_roots_array(data->start, k_perp_samples, data->k_perps, data->omegas, count);
    data->length = N;        // Store the number of calculated roots in the struct

    return NULL;
}


void write_values()
{
    FILE * fout = fopen(VALUEFILE, "w");

    fprintf(fout, "{\n\t\"LAMBDA_C\": %.2f,\n\t\"LAMBDA_H\": %.2f,\n\t\"KAPPA_C\": %.1f,\n\t\"KAPPA_H\": %.1f,\n\t\"N0H_BY_N0E\": %.1f,\n\t\"TH_BY_TC\": %.3f\n}\n", LAMBDA_C, LAMBDA_H, KAPPA_C, KAPPA_H, N0H_BY_N0E, TH_BY_TC);
}
