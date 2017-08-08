# Common utilities used by scripts used to generate series of data for plotting
# by modifying constants.h and running the program repeatedly
#
# We use the sh.py library to issue the system (bash) commands needed to get the desired output

import sh
import sys

SUFFICES = 'abcdef'     # Suffices to be appended to the plot name


# Create backed commands for repeated use
roots = sh.sed.bake('-i', 'roots.h', '-e')
constants = sh.sed.bake('-i', 'constants.h', '-e')        # constants.('foo') will now execute as 'sed -i constants.h -e foo'



def set_defaults(defaults):
    """
    Set default constant values by passing in a dictionary relating variable names to values.
    """

    # Use sed to set the default values
    for k, v in defaults.items():
        constants('s/{key}.*$/{key} {value}/'.format(key=k, value=v))



def set_roots_value(var, value):
    """
    Set the specified value in the roots.h file.
    """

    roots('s/{key}.*$/{key} {value}/'.format(key=var, value=value))



def iterate_variable(plot, var, values):
    """
    Iterate over the values of the specified variable (var) and generate data.

    plot: Name/Number of the plot this data is associated with.
    """

    for i in range(len(values)):

        print("Calculating for {var} = {value} ...\n".format(var=var, value=values[i]))

        constants('s/PLOT.*$/PLOT "{plot}-{suffix}"/'.format(plot=plot, suffix=SUFFICES[i]))
        constants('s/{var}.*$/{var} {value}/'.format(var=var, value=values[i]))

        sh.make('data', _out=sys.stdout, _err=sys.stderr)


        ## Carry out some cleaning on the created data
        #sh.sed('-i', 'data/value-{plot}-{suffix}.json'.format(plot=PLOT, suffix=SUFFICES[i]), '-e', 's/inf/"inf"/')        # Strings need to be quoted in json


    print("Done")
