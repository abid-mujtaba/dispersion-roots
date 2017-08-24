/*
 * Declare all constants used by the numerical recipe.
 */

#include <math.h>
#include <mpfr.h>

// Define the number and version of the plot whose data is being calculated
#define PLOT "00-H"


// Set to 1 if you want debug information printed
#define DEBUG 1

// Define the constants that control the program logic
#define RND MPFR_RNDN           // Set the type of the rounding when using the RND
#define MIN_PRECISION 128            // Bits of precision for MPFR floats
#define TOLERANCE 1e-20                // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function

// If tolderance is NOT achieved the summation will be truncated at this many terms.
// Usually 1000. Needs to be 100,000 for the Kappa -> Infinity case (where K_PERP_MAX = 100)
#define MAX_TERMS 1000




// We use macros to define all variables whose values we would need to adjust. The rest are derived from these.
#define KAPPA_C 2.0
#define KAPPA_H 4.0

// Capital Lambda used in VC Dispersion
#define LAMBDA 0.5

// Relative density of hot electrons out of total (hot and cold) electrons)
#define N0H_BY_N0E 1.0

#define OMEGA_UH_BY_OMEGA_CE 5.099
#define TH_BY_TC 101.695

// In the header file we only declare (NOT define) the variables that are required by other functions.

// Normalization constants
#define OMEGA_CE 1
#define RHO_H 1


// Derived constants
#define OMEGA_CC OMEGA_CE
#define OMEGA_CH OMEGA_CE


// Define data structure for encapsulating constants on a per specie basis
struct Constants {
        double rho;
        double n0_by_n0e;
        double omega_c;

        mpfr_t kappa;

        mpfr_t lambda_vc_p2;
        mpfr_t omega_by_omega_c;
        mpfr_t two_lambda;

        mpfr_t csc;
        mpfr_t pi;
        mpfr_t sqrt_pi;

        mpfr_t g_k_p1_2;                // Gamma(kappa + 1/2)
        mpfr_t g_k_m1_2;
        mpfr_t g_k_p1;
        mpfr_t g_k_m1;
        mpfr_t g_k_p3_2;
        mpfr_t g_k_m3_2;                // Gamma(kappa - 3/2)

        mpfr_t g_2k;                    // Gamma(2 * kappa)
        mpfr_t g_m1_2_mk;               // Gamma(-1/2 - kappa)
        mpfr_t g_1_2_mk;                // Gamma(1/2 - kappa)
        mpfr_t g_5_2_mk;                // Gamma(5/2 - kappa)
};


void get_constants_h(struct Constants * const c);
void get_constants_c(struct Constants * const c);

void clear_constants(struct Constants * const c);
