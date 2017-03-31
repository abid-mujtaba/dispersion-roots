/*
 * Define root finding function prototypes and related constants.
 */

// Define the bracket limits for root finding
#define ROOT_LO 0           // Function is nan at zero so we avoid it
#define ROOT_HI 100
#define ROOT_INTERVAL 1e-3
#define ROOT_MAX_ITERATIONS 13          // to get to 1e-3 from an initial interval of 5 takes about 13 binary divisions

// Define the threshold interval for finding the root. Once the bracket becomes
// smaller than this value the root-finding iteration will stop
#define ROOT_BRACKET 1e-4


// ToDo: Remove prototype (being used in test)
double find_k_perp_root(const double omega, const double lo, const double hi);
int find_k_perp_roots_array(double slices[], double omega[], double roots[], int size);
