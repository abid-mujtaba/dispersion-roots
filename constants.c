/*
 * Populate data structures encapsulating the constants (on a per specie basis)
 */

#include "constants.h"
#include "derived.h"


void get_constants_j(struct Constants * c, double rho, double n0j_by_n0e, double kappa);


void get_constants_h(struct Constants * const c)
{
        get_constants_j(c, RHO_H, N0H_BY_N0E, KAPPA_H);
}


void get_constants_c(struct Constants * const c)
{
        get_constants_j(c, RHO_C, N0C_BY_N0E, KAPPA_C);
}


void get_constants_j(struct Constants * const c, double rho, double n0j_by_n0e, double kappa)
{
        c->rho = rho;
        c->n0_by_n0e = n0j_by_n0e;

        mpfr_init(c->kappa);

        mpfr_set_d(c->kappa, kappa, RND);
}


void clear_constants(struct Constants * const c)
{
        mpfr_clear(c->kappa);
}
