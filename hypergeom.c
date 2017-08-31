/*
 * Define the 1F2 and 2F3 generalized hypergeometric functions required by the dispersion relation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "constants.h"
#include "hypergeom.h"


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

                mpfr_add_d(v, * c.a1, k, RND);
                mpfr_mul(term, term, v, RND);
                mpfr_add_d(v, * c.a2, k, RND);
                mpfr_mul(term, term, v, RND);

                mpfr_mul(term, term, x, RND);

                mpfr_add_d(v, * c.b1, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, * c.b2, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, * c.b3, k, RND);
                mpfr_div(term, term, v, RND);

                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_fprintf(stderr, "\nWarning: Summation truncated before tolerance was achieved. 2F3(%RG, %RG; %RG, %RG, %RG; %RG) - last term = %RG", c.a1, c.a2, c.b1, c.b2, c.b3, x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}


void hyp2F2(mpfr_t result, const struct coeffs_2f2 c, const mpfr_t x)
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

                mpfr_add_d(v, * c.a1, k, RND);
                mpfr_mul(term, term, v, RND);
                mpfr_add_d(v, * c.a2, k, RND);
                mpfr_mul(term, term, v, RND);

                mpfr_mul(term, term, x, RND);

                mpfr_add_d(v, * c.b1, k, RND);
                mpfr_div(term, term, v, RND);
                mpfr_add_d(v, * c.b2, k, RND);
                mpfr_div(term, term, v, RND);

                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_fprintf(stderr, "\nWarning: Summation truncated before tolerance was achieved. 2F2(%RG, %RG; %RG, %RG; %RG) - last term = %RG", c.a1, c.a2, c.b1, c.b2, x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}


// Normalized by dividing by Gamma(c.b3)
void norm_hyp2F3(mpfr_t result, const struct coeffs_2f3 c, const mpfr_t x)
{
        mpfr_t term, fterm, v;
        mpfr_inits(term, fterm, v, (mpfr_ptr) 0);

        mpfr_set_d(result, 0 , RND);
        mpfr_set_d(term, 1, RND);         // Stores the running value of each term in the summation. From the definition of the Pochhammer symbols the value of the k = 0 terms is ONE and then we divide it by the gamma functions

        double b3 = mpfr_get_d(* c.b3, RND);
        int start = 0;

        // We check if c.b3 is possibly a negative integer (or zero)
        if ((b3 <= 0) & (b3 == (int) b3))
                start = (-b3) + 1;              // In this case the first few terms will have to be neglected since they will be divided by the Gamma of a non-positive integer which is infinity and so upon division these terms will become zero


        // We add 'start' to c.b2. If it is zero, no harm. If it is not zero then c.b2 + start = 1 and Gamma(1) = 1 so still no harm
        mpfr_add_ui(v, * c.b3, start, RND);
        mpfr_gamma(v, v, RND);
        mpfr_div(term, term, v, RND);

        mpfr_abs(fterm, term, RND);

        int k;

        for (k = 0; (k < MAX_TERMS) & (mpfr_cmp_d(fterm, TOLERANCE) > 0); ++k)
        {
                // for start != 0 the first 'start' terms are to be neglected so we don't add them to the result
                if (k >= start)
                        mpfr_add(result, result, term, RND);

                mpfr_add_ui(v, * c.a1, k, RND);
                mpfr_mul(term, term, v, RND);

                mpfr_add_ui(v, * c.a2, k, RND);
                mpfr_mul(term, term, v, RND);

                // We use the fact that \Gamma(x + 1) = x * \Gamma(x) to incorporate the increasing gamma values in the term calculation
                mpfr_add_ui(v, * c.b1, k, RND);
                mpfr_div(term, term, v, RND);

                mpfr_add_ui(v, * c.b2, k, RND);
                mpfr_div(term, term, v, RND);

                if (k >= start)         // For start != 0 Gamma(c.b2 + k) = infinity and so the term is neglected
                {
                        mpfr_add_ui(v, * c.b3, k, RND);
                        mpfr_div(term, term, v, RND);
                }

                mpfr_mul(term, term, x, RND);
                mpfr_div_d(term, term, k + 1, RND);

                mpfr_abs(fterm, term, RND);
        }

        if (DEBUG & (k == MAX_TERMS))
            mpfr_fprintf(stderr, "\nWarning: Summation truncated before tolerance was achieved. norm_2F3(x = %RG) - last term = %RG", x, term);

        mpfr_clears(term, fterm, v, (mpfr_ptr) 0);
}
