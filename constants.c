/*
 * Populate data structures encapsulating the constants (on a per specie basis)
 */

#include "constants.h"
#include "derived.h"


void get_constants_j(struct Constants * c, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2);


void get_constants_h(struct Constants * const c)
{
        get_constants_j(c, KAPPA_H, RHO_H, N0H_BY_N0E, OMEGA_CH, LAMBDA_VCH_P2);
}


void get_constants_c(struct Constants * const c)
{
        get_constants_j(c, KAPPA_C, RHO_C, N0C_BY_N0E, OMEGA_CC, LAMBDA_VCC_P2);
}


void get_constants_j(struct Constants * const c, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2)
{
        c->rho = rho;
        c->n0_by_n0e = n0j_by_n0e;
        c->omega_c = omega_c;

        mpfr_inits(c->kappa, c->lambda_vc_p2, c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, (mpfr_t *) 0);

        mpfr_set_d(c->kappa, kappa, RND);
        mpfr_set_str(c->lambda_vc_p2, str_lambda_vc_p2, 10, RND);

        mpfr_const_pi(c->pi, RND);
}


void clear_constants(struct Constants * const c)
{
        mpfr_clears(c->kappa, c->lambda_vc_p2,  c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, (mpfr_t *) 0);
}
