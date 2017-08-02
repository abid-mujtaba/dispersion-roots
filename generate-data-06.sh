#! /bin/bash

# Script that is used to generate series of data for the first plot
# by modifying constants.h and running the program repeatedly

# First we make the global changes (values that are constant across the entire series)
sed -i constants.h  -e 's/KAPPA_C.*$/KAPPA_C 2.0/' \
                    -e 's/KAPPA_H.*$/KAPPA_H 4.0/' \
                    -e 's/N0H_BY_N0E.*$/N0H_BY_N0E 0.0/' \
                    -e 's/LAMBDA.*$/LAMBDA 0.15/'


# Set the value of K_PERP_MAX
sed -i roots.h -e 's/K_PERP_MAX.*$/K_PERP_MAX 100/'



# Now we modify the value of LAMBDA and generate the corresponding data
sed -i constants.h  -e 's/PLOT.*$/PLOT "06-a"/' \
                    -e 's/TH_BY_TC.*$/TH_BY_TC 10.0/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "06-b"/' \
                    -e 's/TH_BY_TC.*$/TH_BY_TC 100.0/'
make data

sed -i constants.h  -e 's/PLOT.*$/PLOT "06-c"/' \
                    -e 's/TH_BY_TC.*$/TH_BY_TC 1000.0/'
make data



# Beep three times to indicate end of the calculation
sudo ~/bin/beeps
