/*
 * Populate data structures encapsulating the constants (on a per specie basis)
 */

#include "constants.h"
#include "derived.h"


void get_constants_j(struct Constants * c, double lambda, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2);


void get_constants_h(struct Constants * const c)
{
        get_constants_j(c, LAMBDA_H, KAPPA_H, RHO_H, N0H_BY_N0E, OMEGA_CH, LAMBDA_VC_P2_H);
}


void get_constants_c(struct Constants * const c)
{
        get_constants_j(c, LAMBDA_C, KAPPA_C, RHO_C, N0C_BY_N0E, OMEGA_CC, LAMBDA_VC_P2_C);
}


void get_constants_j(struct Constants * const c, double lambda, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2)
{
        c->rho = rho;
        c->n0_by_n0e = n0j_by_n0e;
        c->omega_c = omega_c;

        mpfr_inits(c->lambda, c->kappa, c->lambda_vc_p2, c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, c->sqrt_pi, (mpfr_t *) 0);

        mpfr_set_d(c->kappa, kappa, RND);
        mpfr_set_d(c->lambda, lambda, RND);
        mpfr_set_str(c->lambda_vc_p2, str_lambda_vc_p2, 10, RND);

        mpfr_const_pi(c->pi, RND);
        mpfr_sqrt(c->sqrt_pi, c->pi, RND);
}


void clear_constants(struct Constants * const c)
{
        mpfr_clears(c->lambda, c->kappa, c->lambda_vc_p2,  c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, c->sqrt_pi, (mpfr_t *) 0);
}
