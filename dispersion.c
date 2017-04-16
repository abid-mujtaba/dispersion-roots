/*
 * Define the various mathematical functions that will be required to construct
 * the dispersion relation whose roots need to be found.
 */

#include <stdio.h>
#include <mpfr.h>
#include <gsl/gsl_sf_gamma.h>
#include "dispersion.h"
#include "henning.h"            // For lambda_kappa_j_p2(,,)
#include "constants.h"


// Function prototypes
double specie_j(double k_perp, double omega, double lambda_kappa_j_p2, double kappa_j, double omega_cj, double rho_j);
// double specie_j_zero(double kappa_j, double rho_j, struct coeffs_2f3 c_2f3, double lambda_kappa_j_p2);
double specie_c(double k_perp, double omega);
double specie_h(double k_perp, double omega);

double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e);
void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp);


// The following functions are declared here but are defined elsewhere
void calc_first(mpfr_t first, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi);
void calc_second(mpfr_t second, mpfr_t kappa, mpfr_t omega_by_omega_cj, mpfr_t two_lambda_j, mpfr_t csc, mpfr_t pi);


double D(const double k_perp, const double omega)
{
        int p = 0;

        // We start by setting the default precision for MPFR variables based on the value of k_perp. The larger it is the higher the precision required.
        p = 1 + (int) (k_perp / 30);

        mpfr_set_default_prec(MIN_PRECISION * (int) pow(2, p));

        double r = 1 + (specie_c(k_perp, omega) + specie_h(k_perp, omega));

        mpfr_free_cache();              // Needs to be called when constants (like pi have been calculated)


        // Testing using printf statements
        // printf("\nD(%.2f, %.2f) = %.17g", k_perp, omega, r);


        return r;
}


double specie_c(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, KAPPA_C, OMEGA_CC, RHO_C, N0C_BY_N0E);
}

double specie_h(const double k_perp, const double omega)
{
        return specie_j(k_perp, omega, KAPPA_H, OMEGA_CH, RHO_H, N0H_BY_N0E);
}


double specie_j(const double k_perp, const double omega, const double kappa_j, const double omega_cj, const double rho_j, const double n0j_by_n0e)
{
        // ToDo: Deal with the special case k_perp == 0

        double r;

        mpfr_t result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, first, second, third;
        mpfr_inits(result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, first, second, third, (mpfr_ptr) 0);

        // Calculate MPFR variables required for the three terms
        mpfr_set_d(kappa, kappa_j, RND);
        calc_omega_by_omega_cj(omega_by_omega_cj, omega, omega_cj);
        calc_two_lambda_j(two_lambda_j, kappa, rho_j, k_perp);
        mpfr_const_pi(pi, RND);

        mpfr_mul(csc, omega_by_omega_cj, pi, RND);    // csc = omega_by_omega_cj * pi
        mpfr_csc(csc, csc, RND);                      // csc = cosec( csc )


        calc_first(first, kappa, omega_by_omega_cj, two_lambda_j, csc, pi);
        calc_second(second, kappa, omega_by_omega_cj, two_lambda_j, csc, pi);
        // ToDo: Remove this place-holder initial value
        mpfr_set_d(third, 0, RND);


        mpfr_sub(result, first, second, RND);         // result = first - second
        mpfr_sub(result, result, third, RND);         // result -= third

        // Final division
        mpfr_div_d(result, result, pow(k_perp, 2) * lambda_vcj_p2(kappa_j, rho_j, n0j_by_n0e), RND);         // result /= pow(k_perp, 2) * lambda_vcj_p2

        r = mpfr_get_d(result, RND);

        mpfr_clears(result, kappa, omega_by_omega_cj, two_lambda_j, pi, csc, first, second, third, (mpfr_ptr) 0);
        mpfr_free_cache();              // Clear the creation of the constant pi

        return r;
}

/*
 * Calculate specie_j for the special case of k_perp = 0.
 * In this case the only term that survives is the first term of 2F3, all other terms going to zero.
 */
// double specie_j_zero(double kappa_j, double rho_j, struct coeffs_2f3 c, double lambda_kappa_j_p2)
// {
//         // ToDo: Implement
//         return 0;
// }



void calc_two_lambda_j(mpfr_t result, const mpfr_t kappa_j, const double rho_j, const double k_perp)
{
        mpfr_t x;
        mpfr_init(x);


        mpfr_set_d(result, k_perp, RND);                // result = k_perp;
        mpfr_mul_d(result, result, rho_j, RND);         // result *= rho_j;
        mpfr_pow_ui(result, result, 2, RND);            // result = pow(result, 2);

        mpfr_sub_d(x, kappa_j, 1.5, RND);               // x = kappa - 1.5
        mpfr_mul_ui(x, x, 2, RND);                      // x *= 2

        mpfr_mul(result, result, x, RND);               // result *= x


        mpfr_clear(x);
}


double lambda_vcj_p2(double kappa_j, double rho_j, double n0j_by_n0e)
{
        return (3 * LAMBDA + 1) * (kappa_j + 1) * gsl_sf_gamma(kappa_j + 1.5) * lambda_kappa_j_p2(kappa_j, rho_j, n0j_by_n0e);
}
