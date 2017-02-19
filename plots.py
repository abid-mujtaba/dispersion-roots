"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt
import plac

c_functions = ctypes.CDLL('./libfunctions.so')          # Load the library

c_gamma_n_array = c_functions.Gamma_n_array
c_I_n_array = c_functions.I_n_array
c_D_array = c_functions.D_array

c_gamma_n_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_gamma_n_array.restype = None

c_I_n_array.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_I_n_array.restype = None

c_D_array.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_D_array.restype = None


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


def plot_Gamma_n():
    """
    Plot graphs of Gamma_1, Gamma_2, and Gamma_3 from 0 to 4 on the same plot.
    """

    xs = [x / 100.0 for x in range(4 * 100)]

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

    xs = [x / 100.0 for x in range(int(3.5 * 100))]

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


def plot_D():
    """
    Plot a graph of D(omega) (the dispersion relation).
    """

    xs = [x / 1000.0 for x in range(int(4.75 * 1000))]

    Ds = D_array(xs)

    plt.plot(xs, Ds)

    plt.xlabel("$\omega$")
    plt.ylabel("$\mathbf{D}(k_\perp, \omega)$", rotation=0)

    axes = plt.gca()
    axes.set_ylim([-10, 10])


def main(gamma_n: ("Plot Gamma_n for n = 1,2,3", "flag", "g"),
         I_n: ("Plot I_n for n = 0,1,2", "flag", "i")):

    if gamma_n:
        plot_Gamma_n()

    elif I_n:
        plot_I_n()

    # Default option
    else:
        # plot_Gamma_n()
        # plot_I_n()
        plot_D()

    plt.legend()
    plt.grid(True)

    plt.show()


if __name__ == '__main__':
    plac.call(main)
