#! /bin/bash

# Script that is used to generate series of data for the first plot
# by modifying constants.h and running the program repeatedly

# First we make the global changes (values that are constant across the entire series)
sed -i constants.h  -e 's/KAPPA_C.*$/KAPPA_C 2.0/' \
                    -e 's/KAPPA_H.*$/KAPPA_H 4.0/' \
                    -e 's/N0H_BY_N0E.*$/N0H_BY_N0E 0.0/'

# Set the value of K_PERP_MAX
sed -i roots.h -e 's/K_PERP_MAX.*$/K_PERP_MAX 100/'


# Now we modify the value of LAMBDA and generate the corresponding data
sed -i constants.h  -e 's/PLOT.*$/PLOT "01-b-a"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.00/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-b-b"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.10/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-b-c"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.20/'
make data


# Beep three times to indicate end of the calculation
sudo ~/bin/beeps
