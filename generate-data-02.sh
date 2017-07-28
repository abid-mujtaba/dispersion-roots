#! /bin/bash

# Script that is used to generate series of data for the first plot
# by modifying constants.h and running the program repeatedly

# First we make the global changes (values that are constant across the entire series)
sed -i constants.h  -e 's/LAMBDA.*$/LAMBDA 0.15/'  \
                    -e 's/KAPPA_C.*$/KAPPA_C 2.0/' \
                    -e 's/N0H_BY_N0E.*$/N0H_BY_N0E 0.5/'

sed -i roots.h -e 's/K_PERP_MAX.*$/K_PERP_MAX 20/'


# Now we modify the value of KAPPA_H and generate the corresponding data
sed -i constants.h  -e 's/PLOT.*$/PLOT "02-a"/' \
                    -e 's/KAPPA_H.*$/KAPPA_H 1.6/'
make data


sed -i constants.h  -e 's/PLOT.*$/PLOT "02-b"/' \
                    -e 's/KAPPA_H.*$/KAPPA_H 2.0/'
make data


sed -i constants.h  -e 's/PLOT.*$/PLOT "02-c"/' \
                    -e 's/KAPPA_H.*$/KAPPA_H INFINITY/'
make data


# Carry out some cleaning on the created data
sed -i data/value-02-c.json -e 's/inf/"inf"/'       # strings need to be quoted in json


# Beep three times to indicate end of the calculation
sudo ~/bin/beeps
