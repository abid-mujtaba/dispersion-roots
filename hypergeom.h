#define TOLERANCE 1e-20                // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function
#define MAX_TERMS 150                   // If tolderance is NOT achieved the summation will be truncated at this many terms

double hyp1F2(double a1, double b1, double b2, double x);
double hyp2F3(double a1, double a2, double b1, double b2, double b3, double x);

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

double series_hyp(const double coeff, const struct coeffs_1f2 c_1f2, const struct coeffs_2f3 c_2f3, const double x);
double together_hyp(const double coeff, const struct coeffs_1f2 c_1f2, const struct coeffs_2f3 c_2f3, const double x);
