#! /usr/bin/env python3
#
# Script for generating data for plot 07 (different values of KAPPA_H)

# Investigate the dependence of cross-over point on kappa in Henning (single-species)


from generate_data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "07"
VARIABLES = {
            'KAPPA_H': [1.6, 4.0, 10.0],
    }

# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA': 0.00,
        'KAPPA_C': 2.0,
        'N0H_BY_N0E': 1.0,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 30



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
