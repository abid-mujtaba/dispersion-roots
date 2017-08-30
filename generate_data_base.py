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
data = sh.sed.bake('-i', 'data.c', '-e')



def substitute(cmd, var, value):
    """
    cmd: One of the baked commands: roots, constants or data
    var: Variable to focus on
    value: Value to set for spcified variable
    """

    # Note: We include #define so that macro substitution still works (especially in data.c)
    cmd('s/#define {key}.*$/#define {key} {value}/'.format(key=var, value=value))



def set_defaults(defaults):
    """
    Set default constant values by passing in a dictionary relating variable names to values.
    """

    # Use sed to set the default values
    for k, v in defaults.items():
        substitute(constants, k, v)



def set_roots_value(var, value):
    """
    Set the specified value in the roots.h file.
    """

    substitute(roots, var, value)



def set_data_value(var, value):
    """
    Set the specified value in the data.c file.
    """

    substitute(data, var, value)



def iterate_variables(plot, vars):
    """
    Iterate over the values of the specified variables (vars dict) and generate data.

    plot: Name/Number of the plot this data is associated with.
    """

    # We access (arbitrarily) the first element of the dict to get the length of the corresponding dictionary
    L = len(vars[list(vars.keys())[0]])

    for i in range(L):
        for var in vars.keys():

            print("Calculating for {var} = {value} ...".format(var=var, value=vars[var][i]))

            constants('s/{var}.*$/{var} {value}/'.format(var=var, value=vars[var][i]))

        print()
        constants('s/PLOT.*$/PLOT "{plot}-{suffix}"/'.format(plot=plot, suffix=SUFFICES[i]))

        sh.make('data', _out=sys.stdout, _err=sys.stderr)


        ## Carry out some cleaning on the created data
        #sh.sed('-i', 'data/value-{plot}-{suffix}.json'.format(plot=PLOT, suffix=SUFFICES[i]), '-e', 's/inf/"inf"/')        # Strings need to be quoted in json


    print("Done")
