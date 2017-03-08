#ifndef HYPERGEOM_H
#define HYPERGEOM_H


#define TOLERANCE 1e-20                 // Tolerance to be achieved by successive values of the sum while calculating the hypergeometric function
#define MAX_TERMS 500                   // If tolderance is NOT achieved the summation will be truncated at this many terms

double hyp1F2(double a1, double b1, double b2, double x);
double hyp2F3(double a1, double a2, double b1, double b2, double b3, double x);


#endif
