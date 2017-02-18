"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library

_cfunctions.In_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.In_array.restype = None


def In_array(n, xs):

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_Ins = array_type()    # Create a new C-type array that has the capacity to carry the result back
    # Note: We didn't specify any list as the argument so an empty C-type array is created

    _cfunctions.In_array(n, array_type(*xs), c_Ins, num)        # Note how *xs needs to be cast to ctype but n and num can be sent as is and is automatically cast

    return c_Ins            # Python seems capable of handling the c-type array without converting to a standard Python list


def main():

    xs = [x / 100.0 for x in range(500)]
    I0s = In_array(0, xs)

    plt.plot(xs, I0s)
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    main()
