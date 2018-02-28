#! /usr/bin/env python
#
# Script for generating data for plot 10 (compare Maxwellian, Double-Kappa, Cairns and VC dispersion for cold electrons)


from data_base import set_defaults, set_roots_value, iterate_variables


PLOT = "10"
VARIABLES = {
	'LAMBDA_C': [0.0, 0.0, 0.05, 0.05],
	'LAMBDA_H': [0.0, 0.0, 0.2, 0.2],
	'KAPPA_C': ['INFINITY', 2.0, 'INFINITY', 2.0],
	'KAPPA_H': ['INFINITY', 4.0, 'INFINITY', 4.0],
}


# We first declare the default values we want to set as a dictionary
defaults = {
        'N0H_BY_N0E': 0.0,
        'TH_BY_TC': 101.695,
    }

K_PERP_MAX = 100



# Apply changes and generate data
set_roots_value("K_PERP_MAX", K_PERP_MAX)
set_defaults(defaults)
iterate_variables(PLOT, VARIABLES)
