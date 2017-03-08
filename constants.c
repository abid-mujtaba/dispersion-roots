/*
 * Define the constants used externally.
 */

#include "constants.h"
#include <math.h>

// Normalization constants
#define CONST_OMEGA_CE 1
#define CONST_RHO_H 1


// Define physical constants (from NIST database)
const double EPSILON_0 = 8.854187817e-12;         // SI Units
const double ELECTRON_CHARGE = 1.6021766208-19;   // Coulombs
const double ELECTRON_MASS = 9.10938356-31;       // Kilograms

// Derived constants
const double OMEGA_CC = CONST_OMEGA_CE;
const double OMEGA_CH = CONST_OMEGA_CE;

/*
 *Since we are declaring RHO_C as 'const' it can only be constructed from constant objects (not const objects) and so we use a CONST_RHO_H macro to ensure that the same value is used to calculate both RHO_H and RHO_C
*/
const double RHO_H = CONST_RHO_H;
const double RHO_C = CONST_RHO_H / sqrt(TH_BY_TC);
