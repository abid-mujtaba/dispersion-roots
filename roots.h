/*
 * Define root finding function prototypes and related constants.
 */

#ifndef ROOTS_H
#define ROOTS_H

// Define the bracket limits for root finding
#define ROOT_LO 0.0
#define ROOT_HI 5.5
#define ROOT_INTERVAL 1e-3
#define ROOT_MAX_ITERATIONS 13          // to get to 1e-3 from an initial interval of 5 takes about 13 binary divisions

// Define the threshold interval for finding the root. Once the bracket becomes
// smaller than this value the root-finding iteration will stop
#define ROOT_BRACKET 1e-4


int find_k_perp_roots_array(double slices[], double omega[], double roots[], int size);

#endif
