#! /usr/bin/env python
#
# Script for generating data for plot 13 (different values of LAMBDA_H)


from generate_data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "13"
VARIABLES = {
    'LAMBDA_H': [0.00, 0.01, 0.10, 0.25],
}


# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA_C': 0.01,
        'KAPPA_C': 2.0,
        'KAPPA_H': 4.0,
        'N0H_BY_N0E': 0.5,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 10



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
