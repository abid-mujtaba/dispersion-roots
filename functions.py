"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
# from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy


c_functions = ctypes.CDLL('./libDroots.so')          # Load the library

c_D = c_functions.D
c_D_roots = c_functions.find_k_perp_roots_array

c_D.argtypes = (ctypes.c_double, ctypes.c_double)
c_D.restype = ctypes.c_double

c_D_roots.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_D_roots.restype = ctypes.c_int


def D_mesh(k_perp, omega):
    """
    The function takes two arrays of values for k_perp and omega that define the
    coordinates in the mesh.

    It returns a 2D numpy array with the function c_D applied to these coordinates on the mesh.
    """

    D = numpy.zeros(shape=(len(omega), len(k_perp)))        # Create empty 2D numpy array

    for i in range(len(omega)):
        for j in range(len(k_perp)):

            D[i][j] = c_D(k_perp[j], omega[i])      # Use c_D to populate the array

    return D


def D_roots(slices):
    """
    Caclulate the k_perp root, if present, for the specified values of omega (in slices).

    Return both the omega and k_perp values where roots are found.
    """

    num = len(slices)
    array_type = ctypes.c_double * (num * 2)      # there can be a max of 2 roots per slice so the returning array needs to be larger

    c_ks = array_type()
    c_os = array_type()

    size = c_D_roots(array_type(*slices), c_os, c_ks, num)

    ks = []
    os = []

    # We have to limit the returned arrays to the size returned by c_D_roots
    for i in range(size):
        ks.append(c_ks[i])
        os.append(c_os[i])

    return ks, os
