#define PRECISION 53            // Bits of precision for MPFR floats
#define TOLERANCE 1e-20                // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function
#define MAX_TERMS 150                   // If tolderance is NOT achieved the summation will be truncated at this many terms


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
