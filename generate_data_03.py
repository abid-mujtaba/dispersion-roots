#! /usr/bin/env python
#
# Script for generating data for plot 03 (different values of KAPPA_H)


from generate_data_base import set_defaults, set_roots_value, iterate_variable


PLOT = "03"
VARIABLE = 'KAPPA_C'
VALUES = [1.6, 2.0, "INFINITY"]

# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA': 0.0,
        'KAPPA_H': 4.0,
        'N0H_BY_N0E': 0.5,
        'TH_BY_TC': 101.695,
        'MAX_TERMS': 100000
    }

K_PERP_MAX = 100



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variable(PLOT, VARIABLE, VALUES)
