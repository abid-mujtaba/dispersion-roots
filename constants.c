/*
 * Populate data structures encapsulating the constants (on a per specie basis)
 */

#include "constants.h"
#include "derived.h"


void get_constants_h(struct Constants * const c)
{
        c->rho = RHO_H;
        c->n0_by_n0e = N0H_BY_N0E;
}


void get_constants_c(struct Constants * const c)
{
        c->rho = RHO_C;
        c->n0_by_n0e = N0C_BY_N0E;
}
