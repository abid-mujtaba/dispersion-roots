"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt
import numpy
import plac

c_functions = ctypes.CDLL('./libfunctions.so')          # Load the library

c_gamma_n_array = c_functions.Gamma_n_array
c_I_n_array = c_functions.I_n_array
c_D_array = c_functions.D_array
c_D = c_functions.D

c_gamma_n_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_gamma_n_array.restype = None

c_I_n_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_I_n_array.restype = None

c_D_array.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_D_array.restype = None

c_D.argtypes = (ctypes.c_double, ctypes.c_double)
c_D.restype = ctypes.c_double


def function_n_array(fn, n, xs):
    """
    fn: A C-function that takes an int 'n', an array of input values and returns
    (by reference) an array of the resulting values.

    This function takes such a function and acts as a wrapper for it.
    """

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_ys = array_type()    # Create a new C-type array that has the capacity to carry the result back
    # Note: We didn't specify any list as the argument so an empty C-type array is created

    fn(n, array_type(*xs), c_ys, num)        # Note how *xs needs to be cast to ctype but n and num can be sent as is and is automatically cast

    return c_ys            # Python seems capable of handling the c-type array without converting to a standard Python list


def Gamma_n_array(n, xs):

    # Use function_n_array to calculate the values by using the C-function
    # 'c_gamma_n_array'
    return function_n_array(c_gamma_n_array, n, xs)


def I_n_array(n, xs):
    return function_n_array(c_I_n_array, n, xs)


def D_array(xs):

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_ys = array_type()    # Create a new C-type array that has the capacity to carry the result back
    # Note: We didn't specify any list as the argument so an empty C-type array is created

    c_D_array(array_type(*xs), c_ys, num)        # Note how *xs needs to be cast to ctype but n and num can be sent as is and is automatically cast

    return c_ys            # Python seems capable of handling the c-type array without converting to a standard Python list


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


def plot_Gamma_n():
    """
    Plot graphs of Gamma_1, Gamma_2, and Gamma_3 from 0 to 4 on the same plot.
    """

    xs = numpy.linspace(0, 4, 4 * 100)

    G1s = Gamma_n_array(1, xs)
    G2s = Gamma_n_array(2, xs)
    G3s = Gamma_n_array(3, xs)

    plt.plot(xs, G1s, 'b-', label="$n = 1$")
    plt.plot(xs, G2s, 'g-', label="$n = 2$")
    plt.plot(xs, G3s, 'k-', label="$n = 3$")

    plt.xlabel(r"$\beta_c = k_\perp \rho_c$")
    plt.ylabel(r"$\Gamma_n$", rotation=0)
    plt.title(r"$\Gamma_n(\beta_c)$")


def plot_I_n():
    """
    Plot graphs of I_0, I_1, and I_2 from 0 to 3.5 on the same plot.
    """

    xs = numpy.linspace(0, 3.5, 3.5 * 100)

    I0s = I_n_array(0, xs)
    I1s = I_n_array(1, xs)
    I2s = I_n_array(2, xs)

    plt.plot(xs, I0s, 'b-', label="$I_0$")
    plt.plot(xs, I1s, 'g-', label="$I_1$")
    plt.plot(xs, I2s, 'k-', label="$I_2$")

    plt.xlabel("$x$")
    plt.ylabel("$I_n(x)$", rotation=0)

    axes = plt.gca()
    axes.set_ylim([-0.25, 3.25])


def plot_D_omega():
    """
    Plot a graph of D(omega) (the dispersion relation).
    """

    xs = numpy.linspace(0, 4.75, 4.75 * 1000)

    Ds = D_array(xs)

    plt.plot(xs, Ds)

    plt.xlabel(r"$\omega$")
    plt.ylabel(r"$\mathbf{D}(k_\perp, \omega)$", rotation=0)

    axes = plt.gca()
    axes.set_ylim([-10, 10])


def plot_D():
    """
    Plot a 2D graph of D(k_perp, omega) (the dispersion relation)
    """

    # Use numpy to create a 1D array of 200 values between 0 and 5
    k_perp = numpy.linspace(0, 5, 200)
    omega = numpy.linspace(0, 5, 200)

    # Use meshgrid to create 2D arrays by repeating the vectors (1D arrays)
    # k_perp and omega over and over again
    # Every element of the 2D K now contains the corresponding value of k_perp
    K, O = numpy.meshgrid(k_perp, omega)      # Create a mesh by performing the outer-product of the two arrays to get a 2D set of tuples
    D = D_mesh(k_perp, omega)      # Apply c_D to the mesh points

    fig, ax = plt.subplots()
    p = ax.pcolor(K, O, D, cmap=plt.cm.RdBu, vmin=abs(D).min(), vmax=abs(D).max())

    fig.colorbar(p)



def main(gamma_n: ("Plot Gamma_n for n = 1,2,3", "flag", "g"),
         I_n: ("Plot I_n for n = 0,1,2", "flag", "i"),
         D_omega: ("Plot D(omega) with k_perp = 1 fixed", "flag", "o")):

    if gamma_n:
        plot_Gamma_n()

    elif I_n:
        plot_I_n()

    elif D_omega:
        plot_D_omega()

    # Default option
    else:
        # plot_Gamma_n()
        # plot_I_n()
        # plot_D_omega()
        plot_D()

    plt.legend()
    plt.grid(True)

    plt.show()


if __name__ == '__main__':
    plac.call(main)
