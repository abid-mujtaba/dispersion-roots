/*
 * Define the constants used externally.
 */

#include "constants.h"
#include <math.h>

// Function prototypes
double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_omega_ce_p2);


// Normalization constants
#define CONST_OMEGA_CE 1
#define CONST_RHO_H 1


// Define physical constants (from NIST database)
#define EPSILON_0 8.854187817e-12         // SI Units
#define ELECTRON_CHARGE 1.6021766208-19   // Coulombs
#define ELECTRON_MASS 9.10938356-31       // Kilograms

// Derived constants
const double OMEGA_CC = CONST_OMEGA_CE;
const double OMEGA_CH = CONST_OMEGA_CE;

/*
 *Since we are declaring RHO_C as 'const' it can only be constructed from constant objects (not const objects) and so we use a CONST_RHO_H macro to ensure that the same value is used to calculate both RHO_H and RHO_C
*/
const double RHO_H = CONST_RHO_H;
const double RHO_C = CONST_RHO_H / sqrt(TH_BY_TC);

#define N0E_BY_OMEGA_CE_p2 (EPSILON_0 * ELECTRON_MASS * (pow(OMEGA_UH_BY_OMEGA_CE, 2) - 1) / pow(ELECTRON_CHARGE, 2))

const double N0C_BY_N0E = 1 - N0H_BY_N0E;       // The sum of the two is 1


// const double LAMBDA_KAPPA_C_p2 = lambda_kappa_j_p2(KAPPA_C, RHO_C, N0C_BY_N0E);
const double LAMBDA_KAPPA_C_p2 = sqrt(EPSILON_0);


double lambda_kappa_j_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return EPSILON_0 * ELECTRON_MASS * (kappa_j - 1.5) * pow(rho_j, 2) / (n0j_by_n0e * N0E_BY_OMEGA_CE_p2 * pow(ELECTRON_CHARGE, 2) * (kappa_j - 0.5));
}
