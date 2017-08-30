#! /usr/bin/env python
#
# Script for generating data for plot 02 (different values of KAPPA_H)


from generate_data_base import set_defaults, set_roots_value, set_data_value, iterate_variables


PLOT = "01"
VARIABLES = {
    'LAMBDA': [0.00, 0.005, 0.01, 0.05, 0.1],
}


# We first declare the default values we want to set as a dictionary
defaults = {
        'KAPPA_C': 2.0,
        'KAPPA_H': 4.0,
        'N0H_BY_N0E': 1.0,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 20
OMEGA_MIN = 3



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_data_value("OMEGA_MIN", OMEGA_MIN)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
