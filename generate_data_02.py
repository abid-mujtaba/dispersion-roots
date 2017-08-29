#! /usr/bin/env python3
#
# Script for generating data for plot 02 (different values of KAPPA_H)


from generate_data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "02"
VARIABLES = {
            'KAPPA_H': [1.6, 2.0, 'INFINITY'],
            'MAX_TERMS': [1000, 1000, 100000],      # The case of kappa_h = inf requires more terms and higher precision
            'MIN_PRECISION': [128, 128, 256],
            'DOUBLE_PRECISION_DELTA': [30, 30, 7],
            'DEBUG': [0, 0, 1],
    }

# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA': 0.15,
        'KAPPA_C': 2.0,
        'N0H_BY_N0E': 0.5,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 30



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
