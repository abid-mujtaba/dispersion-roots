/*
 * Define the prototypes of functions defined in functions.c
 */

// Declare constants

// Omega is being normalized in multiples of omega_c
#define OMEGA_C 1

// The book chooses this ratio of omega_p to omega_c arbitrarily
#define OMEGA_P 2.5 * OMEGA_C

// For numerical plotting we need to choose these values arbitrarily. 1 is always a good choice
#define K_PERP 1
#define RHO_C 1
#define BETA_C K_PERP * RHO_C

// Define related constants.
#define OMEGA_C_2 OMEGA_C * OMEGA_C
#define OMEGA_P_2 OMEGA_P * OMEGA_P
#define BETA_C_2 BETA_C * BETA_C


/* Define the number of terms on either side of the evaluation point to be used
 * while calculating the dispersion relation
 */
#define HALF_MAX_N 10

// Define the bracket limits for root finding
#define ROOT_LO 0.0
#define ROOT_HI 5.0
#define ROOT_INTERVAL 1e-3
#define ROOT_MAX_ITERATIONS 13          // to get to 1e-3 from an initial interval of 5 takes about 13 binary divisions

// Define the threshold interval for finding the root. Once the bracket becomes
// smaller than this value the root-finding iteration will stop
#define ROOT_BRACKET 1e-4


void I_n_array(int n, double x[], double In[], int size);
void Gamma_n_array(int n, double x[], double Gn[], int size);
double D(double k_perp, double omega);
void D_array(double x[], double Ds[], int size);
int find_k_perp_roots(double slices[], double omega[], double roots[], int size);
double find_k_perp_root(double omega);
