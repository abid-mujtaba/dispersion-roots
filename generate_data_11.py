#! /usr/bin/env python
#
# Script for generating data for plot 02 (different values of KAPPA_H)


from generate_data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "11"
VARIABLES = {
    'TH_BY_TC': [10.0, 100.0, 1000.0],
}


# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA_C': 0.2,
        'LAMBDA_H': 0.2,
        'KAPPA_C': 2.0,
        'KAPPA_H': 4.0,
        'N0H_BY_N0E': 0.0,
    }

K_PERP_MAX = 50



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
