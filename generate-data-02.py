#! /usr/bin/env python

# Script that is used to generate series of data for the first plot
# by modifying constants.h and running the program repeatedly
#
# We use the sh.py library to issue the system (bash) commands needed to get the desired output

import sh

PLOT = "02"
VARIABLE = 'KAPPA_H'
VALUES = [1.6, 2.0, 'INFINITY']

# We first declare the default values we want to set as a dictionary
defaults = {
        'LAMBDA': 0.0,
        'KAPPA_C': 2.0,
        'N0H_BY_N0E': 0.5,
        'TH_BY_TC': 101.695
    }

K_PERP_MAX = 20

SUFFICES = 'abcdef'     # Suffices to be appended to the plot name


# Set k_perp max in roots.h
sh.sed('-i', 'roots.h', '-e' 's/K_PERP_MAX.*$/K_PERP_MAX {}/'.format(K_PERP_MAX))


# Create backed commands for repeated use
constants = sh.sed.bake('-i', 'constants.h', '-e')        # constants.('foo') will now execute as 'sed -i constants.h -e foo'


# Use sed to set the default values
for k, v in defaults.items():
    constants('s/{key}.*$/{key} {value}/'.format(key=k, value=v))


for i in range(len(VALUES)):

    print("Calculating for {var} = {value} ...".format(var=VARIABLE, value=VALUES[i]))

    constants('s/PLOT.*$/PLOT "{plot}-{suffix}"/'.format(plot=PLOT, suffix=SUFFICES[i]))
    constants('s/{var}.*$/{var} {value}/'.format(var=VARIABLE, value=VALUES[i]))

    sh.make('data')


## Carry out some cleaning on the created data
sh.sed('-i', 'data/value-{plot}-c.json'.format(plot=PLOT), '-e', 's/inf/"inf"/')        # Strings need to be quoted in json

print("Done")
