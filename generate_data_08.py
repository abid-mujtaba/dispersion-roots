#! /usr/bin/env python3
#
# Script for generating data for plot 08

# Investigate the dependence of cross-over point on omega in Henning (single-species)


from generate_data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "08"
VARIABLES = {
            'KAPPA_H': [4.0],
    }

# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA_C': 0.00,
        'LAMBDA_H': 0.00,
        'KAPPA_C': 2.0,
        'N0H_BY_N0E': 1.0,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 20



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
