"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library
_cfunctions.J0.argtypes = (ctypes.c_double,)        # Explicitly declare the types of the arguments to the C function
_cfunctions.J0.restype = ctypes.c_double            # Explicity declare the type of the result of the function


def J0(x):

    global _cfunctions          # Not strictly necessary when a local variable with the same name doesn't exist

    ans = _cfunctions.J0(ctypes.c_double(x))        # The python object 'x' is cast to a ctypes.c_double type before passing it

    return float(ans)           # The result is explicitly cast back to python's 'float' type


def main():

    x = 5.0 

    print("J0({}) = {}".format(x, J0(x)))


if __name__ == '__main__':
    main()
