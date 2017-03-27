/*
 * Declare all constants used by the numerical recipe.
 */

 #include <math.h>              // Required by the sqrt() functions define below
 #include <mpfr.h>


// Set to 1 if you want debug information printed
#define DEBUG 1

// Define the constants that control the program logic
#define RND MPFR_RNDN           // Set the type of the rounding when using the RND
#define PRECISION 256            // Bits of precision for MPFR floats
#define TOLERANCE 1e-20                // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function
#define MAX_TERMS 1000                   // If tolderance is NOT achieved the summation will be truncated at this many terms




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
void calc_omega_by_omega_cj(mpfr_t result, double omega, double omega_cj);

// Use the function above to define the following
#define LAMBDA_KAPPA_C_p2 lambda_kappa_j_p2(KAPPA_C, RHO_C, N0C_BY_N0E)
#define LAMBDA_KAPPA_H_p2 lambda_kappa_j_p2(KAPPA_H, RHO_H, N0H_BY_N0E)
