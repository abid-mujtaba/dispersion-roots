#! /bin/bash

# Script that is used to generate series of data for the first plot
# by modifying constants.h and running the program repeatedly

# First we make the global changes (values that are constant across the entire series)
sed -i constants.h  -e 's/KAPPA_C.*$/KAPPA_C 2.0/' \
                    -e 's/KAPPA_H.*$/KAPPA_H 4.0/' \
                    -e 's/N0H_BY_N0E.*$/N0H_BY_N0E 1.0/'

# Set the value of K_PERP_MAX
sed -i roots.h -e 's/K_PERP_MAX.*$/K_PERP_MAX 20/'


# Now we modify the value of LAMBDA and generate the corresponding data
sed -i constants.h  -e 's/PLOT.*$/PLOT "01-1-a"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.05/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-1-b"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.10/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-1-c"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.15/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-1-d"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.20/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-1-e"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.25/'
make data


# Change N0H_BY_N0E = 0.0 and repeat series
sed -i constants.h -e 's/N0H_BY_N0E.*$/N0H_BY_N0E 0.0/'


sed -i constants.h  -e 's/PLOT.*$/PLOT "01-2-a"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.05/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-2-b"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.10/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-2-c"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.15/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-2-d"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.20/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "01-2-e"/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.25/'
make data



# Beep three times to indicate end of the calculation
sudo ~/bin/beeps
