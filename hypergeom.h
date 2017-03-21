#include "constants.h"


#define TOLERANCE 1e-20                // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function
#define MAX_TERMS 150                   // If tolderance is NOT achieved the summation will be truncated at this many terms
#define THRESHOLD 125
#define TAYLOR_STEP_K 5                // Step size for shifting center of Taylor expansion in units of k_perp
// Calculate Taylor Step size in 2 \lambda_j^\prime from TAYLOR_STEP_K
// The Hot variants are used since they result in a larger values than the cold variants
#define TAYLOR_STEP (2 * (KAPPA_H - 1.5) * pow(TAYLOR_STEP_K * RHO_H, 2))

struct coeffs_1f2 {
        double a1;
        double b1;
        double b2;
};


struct coeffs_2f3 {
        double a1;
        double a2;
        double b1;
        double b2;
        double b3;
};


int compare_terms(const void *pa, const void *pb);

double hyp1F2(const struct coeffs_1f2 c_1f2, const double x);
double hyp2F3(const struct coeffs_2f3 c_2f3, const double x);
