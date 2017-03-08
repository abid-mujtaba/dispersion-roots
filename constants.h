/*
 * Declare all constants used by the numerical recipe.
 */

// We use macros to define all variables whose values we would need to adjust. The rest are derived from these.
#define KAPPA_C 2
#define KAPPA_H 4

#define N0H_BY_N0E 1.0    // Relative density of hot electrons out of total (hot and cold) electrons

#define OMEGA_UH_BY_OMEGA_CE 5.099
#define TH_BY_TC 101.695

// In the header file we only declare (NOT define) the variables that are required by other functions.
const double OMEGA_CC, OMEGA_CH;
const double RHO_C, RHO_H;
