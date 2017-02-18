"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library

_cfunctions.J0.argtypes = (ctypes.c_double,)        # Explicitly declare the types of the arguments to the C function
_cfunctions.J0.restype = ctypes.c_double            # Explicity declare the type of the result of the function

_cfunctions.sum.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.sum.restype = ctypes.c_double


def J0(x):

    ans = _cfunctions.J0(ctypes.c_double(x))        # The python object 'x' is cast to a ctypes.c_double type before passing it

    return float(ans)           # The result is explicitly cast back to python's 'float' type


def csum(xs):

    num = len(xs)               # C functions require the length of arrays to work
    array_type = ctypes.c_double * num      # Create a new type for an array of type c_double (with the length specified)

    # Note how we expand the list xs using * before passing it to array_type
    return float(_cfunctions.sum(array_type(*xs), ctypes.c_int(num)))


def main():

    x = 5.0
    print("J0({}) = {}".format(x, J0(x)))

    xs = [1.1, 2.2, 3.3, 4.4]
    print("Sum of {} = {}".format(str(xs), csum(xs)))


if __name__ == '__main__':
    main()
