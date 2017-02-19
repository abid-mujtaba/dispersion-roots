"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt
import plac

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library

_cfunctions.Gamma_n_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.Gamma_n_array.restype = None


def Gamma_n_array(n, xs):

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_Gns = array_type()    # Create a new C-type array that has the capacity to carry the result back
    # Note: We didn't specify any list as the argument so an empty C-type array is created

    _cfunctions.Gamma_n_array(n, array_type(*xs), c_Gns, num)        # Note how *xs needs to be cast to ctype but n and num can be sent as is and is automatically cast

    return c_Gns            # Python seems capable of handling the c-type array without converting to a standard Python list


def plot_Gamma_n():
    """
    Plot graphs of Gamma_1, Gamma_2, and Gamma_3 from 0 to 4 on the same plot.
    """

    xs = [x / 100.0 for x in range(400)]

    G1s = Gamma_n_array(1, xs)
    G2s = Gamma_n_array(2, xs)
    G3s = Gamma_n_array(3, xs)

    plt.plot(xs, G1s, 'b-')
    plt.plot(xs, G2s, 'g-')
    plt.plot(xs, G3s, 'k-')
    plt.grid(True)
    plt.show()


def main(gamma_n: ("Plot Gamma_n for n = 1,2,3", "flag", "g")):

    if gamma_n:
        plot_Gamma_n()


if __name__ == '__main__':
    plac.call(main)
