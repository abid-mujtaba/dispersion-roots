/*
 * Define root finding function prototypes and related constants.
 */

// Define the bracket limits for root finding
#define K_PERP_MAX 100          // Max value of k_perp for root finding
#define ROOT_INTERVAL 1e-3
#define DELTA 1e-10             // Minimum value to move away from thresholds (e.g. k_perp = 0 + DELTA since the dispersion relation is undefined at k_perp = 0)


int find_omega_root(const double k_perp, const double lo, const double hi, double * root, const double guess);
int find_omega_roots_array(const int initial, double k_perp_samples[], double k_perps[], double omegas[], const int size);
