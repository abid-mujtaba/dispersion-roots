"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library

_cfunctions.J0_array.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.J0_array.restype = None


def j0_array(xs):

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_j0s = array_type()    # Create a new C-type array that has the capacity to carry the result back
    # Note: We didn't specify any list as the argument so an empty C-type array is created

    _cfunctions.J0_array(array_type(*xs), c_j0s, num)        # Note how *xs needs to be cast to ctype but num can be sent as is and is automatically cast

    return [float(j) for j in c_j0s]        # We expicitly convert to a python list before returning


def main():

    xs = [x / 100.0 for x in range(5000)]
    j0s = j0_array(xs)

    plt.plot(xs, j0s)
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    main()
