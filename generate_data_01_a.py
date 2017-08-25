#! /usr/bin/env python
#
# Script for generating data for plot 02 (different values of KAPPA_H)


from generate_data_base import set_defaults, set_roots_value, iterate_variable


PLOT = "01-a"
VARIABLE = 'LAMBDA'
VALUES = [0.00, 0.10, 0.20]

# We first declare the default values we want to set as a dictionary
defaults = {
        'KAPPA_C': 2.0,
        'KAPPA_H': 4.0,
        'N0H_BY_N0E': 1.0,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 20



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variable(PLOT, VARIABLE, VALUES)
