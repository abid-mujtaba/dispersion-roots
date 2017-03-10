/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include "functions.h"
#include "constants.h"
#include "hypergeom.h"
#include <gsl/gsl_math.h>       // Defines M_SQRTPI
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_gamma.h>


double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j);
double specie_c(double k_perp, double omega);
double specie_h(double k_perp, double omega);

// Place-holder function
double D(const double k_perp, const double omega)
{
        return 1 + (specie_c(k_perp, omega) + specie_h(k_perp, omega));
}


double specie_c(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, LAMBDA_KAPPA_C_p2, KAPPA_C, OMEGA_CC, RHO_C);
}

double specie_h(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, LAMBDA_KAPPA_H_p2, KAPPA_H, OMEGA_CH, RHO_H);
}


double specie_j(const double k_perp, const double omega, const double lambda_kappa_j_p2, const double kappa_j, const double omega_cj, const double rho_j)
{
        const double two_lambda_j_prime = 2 * (kappa_j - 1.5) * pow(k_perp * rho_j, 2);
        const double omega_by_omega_cj = omega / omega_cj;

        double result = 1;
        result -= hyp2F3(1, 0.5, 0.5 - kappa_j, 1 + omega_by_omega_cj, 1 - omega_by_omega_cj, two_lambda_j_prime);

        double third = M_SQRTPI * omega_by_omega_cj;
        third *= gsl_sf_gamma(kappa_j + 1) * gsl_sf_gamma(0.5 - kappa_j);
        third /= gsl_sf_sin(M_PI * omega_by_omega_cj);
        third /= gsl_sf_gamma(kappa_j + 1.5 + omega_by_omega_cj) * gsl_sf_gamma(kappa_j + 1.5 - omega_by_omega_cj);
        third *= pow(two_lambda_j_prime, kappa_j + 0.5);
        third *= hyp1F2(kappa_j + 1, kappa_j + 1.5 + omega_by_omega_cj, kappa_j + 1.5 - omega_by_omega_cj, two_lambda_j_prime);

        result += third;

        result /= (pow(k_perp, 2) * lambda_kappa_j_p2);

        return result;
}
