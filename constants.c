/*
 * Populate data structures encapsulating the constants (on a per specie basis)
 */

#include "constants.h"
#include "derived.h"


void get_constants_j(struct Constants * c, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2, char * str_g_k_p1_2, char * str_g_k_m1_2, char * str_g_k_p1, char * str_g_k_m1, char * str_g_k_p3_2, char * str_g_k_m3_2, char * str_g_2k, char * str_g_m1_2_mk, char * str_g_1_2_mk, char * str_g_5_2_mk);


void get_constants_h(struct Constants * const c)
{
        get_constants_j(c, KAPPA_H, RHO_H, N0H_BY_N0E, OMEGA_CH, LAMBDA_VC_P2_H, G_K_P1_2_H, G_K_M1_2_H, G_K_P1_H, G_K_M1_H, G_K_P3_2_H, G_K_M3_2_H, G_2K_H, G_M1_2_MK_H, G_1_2_MK_H, G_5_2_MK_H);
}


void get_constants_c(struct Constants * const c)
{
        get_constants_j(c, KAPPA_C, RHO_C, N0C_BY_N0E, OMEGA_CC, LAMBDA_VC_P2_C, G_K_P1_2_C, G_K_M1_2_C, G_K_P1_C, G_K_M1_C, G_K_P3_2_C, G_K_M3_2_C, G_2K_C, G_M1_2_MK_C, G_1_2_MK_C, G_5_2_MK_C);
}


void get_constants_j(struct Constants * const c, double kappa, double rho, double n0j_by_n0e, double omega_c, char * str_lambda_vc_p2, char * str_g_k_p1_2, char * str_g_k_m1_2, char * str_g_k_p1, char * str_g_k_m1, char * str_g_k_p3_2, char * str_g_k_m3_2, char * str_g_2k, char * str_g_m1_2_mk, char * str_g_1_2_mk, char * str_g_5_2_mk)
{
        c->rho = rho;
        c->n0_by_n0e = n0j_by_n0e;
        c->omega_c = omega_c;

        mpfr_inits(c->kappa, c->lambda_vc_p2, c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, c->sqrt_pi, c->g_k_p1_2, c->g_k_m1_2, c->g_k_p1, c->g_k_m1, c->g_k_p3_2, c->g_k_m3_2, c->g_2k, c->g_m1_2_mk, c->g_1_2_mk, c->g_5_2_mk, (mpfr_t *) 0);

        mpfr_set_d(c->kappa, kappa, RND);
        mpfr_set_str(c->lambda_vc_p2, str_lambda_vc_p2, 10, RND);

        mpfr_set_str(c->g_k_p1_2, str_g_k_p1_2, 10, RND);
        mpfr_set_str(c->g_k_m1_2, str_g_k_m1_2, 10, RND);
        mpfr_set_str(c->g_k_p1, str_g_k_p1, 10, RND);
        mpfr_set_str(c->g_k_m1, str_g_k_m1, 10, RND);
        mpfr_set_str(c->g_k_p3_2, str_g_k_p3_2, 10, RND);
        mpfr_set_str(c->g_k_m3_2, str_g_k_m3_2, 10, RND);
        mpfr_set_str(c->g_2k, str_g_2k, 10, RND);
        mpfr_set_str(c->g_m1_2_mk, str_g_m1_2_mk, 10, RND);
        mpfr_set_str(c->g_1_2_mk, str_g_1_2_mk, 10, RND);
        mpfr_set_str(c->g_5_2_mk, str_g_5_2_mk, 10, RND);

        mpfr_const_pi(c->pi, RND);
        mpfr_sqrt(c->sqrt_pi, c->pi, RND);
}


void clear_constants(struct Constants * const c)
{
        mpfr_clears(c->kappa, c->lambda_vc_p2,  c->omega_by_omega_c, c->two_lambda, c->csc, c->pi, c->sqrt_pi, c->g_k_p1_2, c->g_k_m1_2, c->g_k_p1, c->g_k_m1, c->g_k_p3_2, c->g_k_m3_2, c->g_2k, c->g_m1_2_mk, c->g_1_2_mk, c->g_5_2_mk, (mpfr_t *) 0);
}
