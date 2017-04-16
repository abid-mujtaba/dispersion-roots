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
        /*
         * For now only the 'term' and 'result' are stored in arbitrary precision.
         * The interim calculations used to calculate term iteratively are still being carried out in double precision
         * Also note that the terms are not being stored in an array and being sorted in the hope that the arbitrary precision will take care of this problem.
         *
         * fterm simply contains the absolutely value of the term
         */
        mpfr_t term, fterm, v;
        mpfr_inits(term, fterm, v, (mpfr_ptr) 0);

        mpfr_set_d(result, 0, RND);
        mpfr_set_d(term, 1, RND);         // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE
        mpfr_set_d(fterm, 1, RND);

        // When (the absolute value of) 'term' becomes less than TOLERANCE the sum stops changing rapidly and we truncate it
        // mpfr_cmp_d returns -1 when fterm < TOLERANCE. Since that is non-zero we need a further > 0 comparison to get the boolean
        // logic to work out.
        int k;

        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                mpfr_add(result, result, term, RND);

                /*
                 * Based on the definition of 1F2 each term differs from the previous by multiplicative factors that have to do with the Pochhammer symbols, power of x and the factorial in the denominator
                 *
                 * Implement the following in MPFR:
                 *      term *= (c.a1 + k) * x / ((c.b1 + k) * (c.b2 + k) * (k + 1));
                 */
                mpfr_add_d(v, c.a1, k, RND);
                mpfr_mul(term, term, v, RND);

                mpfr_mul(term, term, x, RND);

                mpfr_add_d(v, c.b1, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, c.b2, k, RND);
                mpfr_div(term, term, v, RND);

                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_printf("\nWarning: Summation truncated before tolerance was achieved. 1F2(x = %RG) - last term = %RG", x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}


void hyp2F3(mpfr_t result, const struct coeffs_2f3 c, const mpfr_t x)
{
        mpfr_t term, fterm, v;
        mpfr_inits(term, fterm, v, (mpfr_ptr) 0);

        mpfr_set_d(result, 0, RND);
        mpfr_set_d(term, 1, RND);
        mpfr_set_d(fterm, 1, RND);

        int k;

        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                mpfr_add(result, result, term, RND);

                mpfr_add_d(v, c.a1, k, RND);
                mpfr_mul(term, term, v, RND);
                mpfr_add_d(v, c.a2, k, RND);
                mpfr_mul(term, term, v, RND);

                mpfr_mul(term, term, x, RND);

                mpfr_add_d(v, c.b1, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, c.b2, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, c.b3, k, RND);
                mpfr_div(term, term, v, RND);

                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_printf("\nWarning: Summation truncated before tolerance was achieved. 2F3(x = %RG) - last term = %RG", x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}


/* Calculate 1F2 normalized by dividing ONLY by the Gamma function of the b2 variable. */
void norm_hyp1F2(mpfr_t result, const struct coeffs_1f2 c, const mpfr_t x)
{
        mpfr_t term, fterm, v;
        mpfr_inits(term, fterm, v, (mpfr_ptr) 0);

        mpfr_set_d(result, 0, RND);

        mpfr_set_d(term, 1, RND);         // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE and then we divide it by the gamma functions

        // We check if c.b2 is possibly a negative integer (or zero)
        double b2 = mpfr_get_d(c.b2, RND);
        int start = 0;

        if ((b2 <= 0) & (b2 == (int) b2))
                start = (-b2) + 1;              // In this case the first few terms will have to be neglected since they will be divided by the Gamma of a non-positive integer which is infinity and so upon division these terms will become zero


        // We add 'start' to c.b2. If it is zero, no harm. If it is not zero then c.b2 + start = 1 and Gamma(1) = 1 so still no harm
        mpfr_add_d(v, c.b2, start, RND);
        mpfr_gamma(v, v, RND);
        mpfr_div(term, term, v, RND);

        mpfr_abs(fterm, term, RND);


        int k;

        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                // for start != 0 the first 'start' terms are to be neglected so we don't add them to the result
                if (k >= start)
                        mpfr_add(result, result, term, RND);

                mpfr_add_ui(v, c.a1, k, RND);
                mpfr_mul(term, term, v, RND);

                // We use the fact that \Gamma(x + 1) = x * \Gamma(x) to incorporate the increasing gamma values in the term calculation
                mpfr_add_ui(v, c.b1, k, RND);
                mpfr_div(term, term, v, RND);

                if (k >= start)         // For start != 0 Gamma(c.b2 + k) = infinity and so the term is neglected
                {
                        mpfr_add_ui(v, c.b2, k, RND);
                        mpfr_div(term, term, v, RND);
                }

                mpfr_mul(term, term, x, RND);
                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_printf("\nWarning: Summation truncated before tolerance was achieved. 1F2(x = %RG) - last term = %RG", x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}


void init_coeffs(struct coeffs_1f2 * const c1, struct coeffs_2f3 * const c2)
{
        mpfr_inits(c1->a1, c1->b1, c1->b2, c2->a1, c2->a2, c2->b1, c2->b2, c2->b3, (mpfr_ptr) 0);
}


void clear_coeffs(struct coeffs_1f2 * const c1, struct coeffs_2f3 * const c2)
{
        mpfr_clears(c1->a1, c1->b1, c1->b2, c2->a1, c2->a2, c2->b1, c2->b2, c2->b3, (mpfr_ptr) 0);
}
