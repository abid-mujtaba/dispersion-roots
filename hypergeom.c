/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


// Carries out the calculation of 1F2 using MPFR and stores the result in the mpfr_t variable 'result' that is passed in as the first argument
void hyp1F2(mpfr_t result, const struct coeffs_1f2 c, const mpfr_t x)
{
        int k;

        /*
         * For now only the 'term' and 'result' are stored in arbitrary precision.
         * The interim calculations used to calculate term iteratively are still being carried out in double precision
         * Also note that the terms are not being stored in an array and being sorted in the hope that the arbitrary precision will take care of this problem.
         *
         * fterm simply contains the absolutely value of the term
         */
        mpfr_t term, fterm;
        mpfr_inits(term, fterm, (mpfr_ptr) 0);

        mpfr_set_d(result, 0, RND);
        mpfr_set_d(term, 1, RND);         // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE
        mpfr_set_d(fterm, 1, RND);

        // When (the absolute value of) 'term' becomes less than TOLERANCE the sum stops changing rapidly and we truncate it
        // mpfr_cmp_d returns -1 when fterm < TOLERANCE. Since that is non-zero we need a further > 0 comparison to get the boolean
        // logic to work out.
        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                mpfr_add(result, result, term, RND);

                /*
                 * Based on the definition of 1F2 each term differs from the previous by multiplicative factors that have to do with the Pochhammer symbols, power of x and the factorial in the denominator
                 *
                 * Implement the following in MPFR:
                 *      term *= (c.a1 + k) * x / ((c.b1 + k) * (c.b2 + k) * (k + 1));
                 */
                mpfr_mul_d(term, term, c.a1 + k, RND);
                mpfr_mul(term, term, x, RND);
                mpfr_div_d(term, term, c.b1 + k, RND);
                mpfr_div_d(term, term, c.b2 + k, RND);
                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        mpfr_clears(term, fterm, (mpfr_ptr) 0);
}


void hyp2F3(mpfr_t result, const struct coeffs_2f3 c, const mpfr_t x)
{
        int k;

        mpfr_t term, fterm;
        mpfr_inits(term, fterm, (mpfr_ptr) 0);

        mpfr_set_d(result, 0, RND);
        mpfr_set_d(term, 1, RND);
        mpfr_set_d(fterm, 1, RND);

        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                mpfr_add(result, result, term, RND);

                mpfr_mul_d(term, term, c.a1 + k, RND);
                mpfr_mul_d(term, term, c.a2 + k, RND);
                mpfr_mul(term, term, x, RND);

                mpfr_div_d(term, term, c.b1 + k, RND);
                mpfr_div_d(term, term, c.b2 + k, RND);
                mpfr_div_d(term, term, c.b3 + k, RND);
                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        mpfr_clears(term, fterm, (mpfr_ptr) 0);
}
