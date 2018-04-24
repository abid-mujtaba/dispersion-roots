#! /usr/bin/env python
#
# Script for generating data for plot 10 (compare Maxwellian, Double-Kappa, Cairns and VC dispersion for an equal mixture of hot and cold electrons)


from data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "11"
VARIABLES = {
	'LAMBDA_C': [0.0, 0.0, 0.05, 0.05],
	'LAMBDA_H': [0.0, 0.0, 0.2, 0.2],
	'KAPPA_C': ['INFINITY', 3.0, 'INFINITY', 3.0],
	'KAPPA_H': ['INFINITY', 4.0, 'INFINITY', 4.0],
	'MAX_TERMS': [100000, 1000, 100000, 1000],
	'MIN_PRECISION': [256, 128, 256, 128],
	'DOUBLE_PRECISION_DELTA': [7, 30, 7, 30],
}


# We first declare the default values we want to set as a dictionary
defaults = {
        'N0H_BY_N0E': 0.5,
        'TH_BY_TC': 101.695,
    }

K_PERP_MAX = 30



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
