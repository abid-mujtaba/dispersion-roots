/*
 * Declare all constants used by the numerical recipe.
 */

 #include <math.h>              // Required by the sqrt() functions define below

// We use macros to define all variables whose values we would need to adjust. The rest are derived from these.
#define KAPPA_C 2
#define KAPPA_H 4

#define N0H_BY_N0E 1.0    // Relative density of hot electrons out of total (hot and cold) electrons)

#define OMEGA_UH_BY_OMEGA_CE 5.099
#define TH_BY_TC 101.695

// In the header file we only declare (NOT define) the variables that are required by other functions.

// Normalization constants
#define OMEGA_CE 1
#define RHO_H 1


// Derived constants
#define OMEGA_CC OMEGA_CE
#define OMEGA_CH OMEGA_CE

#define RHO_C (RHO_H / sqrt(TH_BY_TC))
#define N0C_BY_N0E (1 - N0H_BY_N0E)       // The sum of the two is 1


// Function prototype
double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_by_n0e);

// Use the function above to define the following
#define LAMBDA_KAPPA_C_p2 lambda_kappa_j_p2(KAPPA_C, RHO_C, N0C_BY_N0E)
#define LAMBDA_KAPPA_H_p2 lambda_kappa_j_p2(KAPPA_H, RHO_H, N0H_BY_N0E)
