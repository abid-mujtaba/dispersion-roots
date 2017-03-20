"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpmath import hyp2f3
import numpy
import plac

c_functions = ctypes.CDLL('./libDroots.so')          # Load the library

# c_D_array = c_functions.D_array
c_D = c_functions.D
c_D_roots = c_functions.find_k_perp_roots_array

# c_D_array.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
# c_D_array.restype = None

c_D.argtypes = (ctypes.c_double, ctypes.c_double)
c_D.restype = ctypes.c_double

c_D_roots.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
c_D_roots.restype = ctypes.c_int


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


def D_roots(slices):
    """
    Caclulate the k_perp root, if present, for the specified values of omega (in slices).

    Return both the omega and k_perp values where roots are found.
    """

    num = len(slices)
    array_type = ctypes.c_double * num

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


def plot_c_array(c_array_function, start=0, end=4, samples=100):
    """
    Generic function for plotting a C function with prototype:
        void c_array_function(double x[], double y[], int size)
    """

    xs = numpy.linspace(start, end, end * samples, endpoint=False)
    num = len(xs)

    array_type = ctypes.c_double * num

    c_xs = array_type(*xs)
    c_ys = array_type()

    c_array_function(c_xs, c_ys, num)

    plt.plot(c_xs, c_ys)


def plot_D_omega():
    """
    Plot a graph of D(omega) (the dispersion relation).
    """

    # Note: Setting endpoint=False means we are NOT including the end-point
    # which will make the intervals even (that is the exact data-points 1,2,3,4
    # will be part of the data set for which the function value is infinity)
    xs = numpy.linspace(0, 4.75, 4.75 * 1000, endpoint=False)

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

    # Use numpy to create a 1D array of 200 values bet8een 0 and 5
    # k_perp = numpy.linspace(0, 100, 1000)
    # omega = numpy.linspace(0, 8, 800)

    k_perp = numpy.linspace(0, 20, 100)
    omega = numpy.linspace(0, 8, 200)

    # Use meshgrid to create 2D arrays by repeating the vectors (1D arrays)
    # k_perp and omega over and over again
    # Every element of the 2D K now contains the corresponding value of k_perp
    K, O = numpy.meshgrid(k_perp, omega)      # Create a mesh by performing the outer-product of the two arrays to get a 2D set of tuples
    D = D_mesh(k_perp, omega)      # Apply c_D to the mesh points

    fig, ax = plt.subplots()

    p = ax.pcolor(K, O, D, cmap=plt.cm.RdBu, vmin=abs(D).min(), vmax=abs(D).max())
    fig.colorbar(p)

    # cnt = plt.contour(K, O, D, cmap=plt.cm.RdBu)
    # fig.colorbar(cnt)


    # 3D mesh surface:

    # fig = plt.figure(figsize=(14, 6))
    #
    # ax = fig.add_subplot(1, 2, 1, projection='3d')
    # p = ax.plot_surface(K, O, D, rstride=4, cstride=4, linewidth=0)
    #
    # ax = fig.add_subplot(1, 2, 2, projection='3d')
    # p = ax.plot_surface(K, O, D, rstride=1, cstride=1, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
    #
    # cb = fig.colorbar(p, shrink=0.5)


def plot_D_roots():
    """
    Plot a graph of the roots of D on axes of omega vs. k_perp.
    """

    for i in range(5):

        slices = numpy.linspace(i, i + 1, 10.0, endpoint=False)

        # We are interested in getting the root very near the boundary but not on it since the function blows up at the boundary.
        # So we append a value to slices which is extremely close to the boundary
        slices = numpy.append(slices, i + 1e-10)

        K, O = D_roots(slices)

        # The data returned is for slices of Omega with up to two posible values of k_perp
        # To get a smooth curve we will have to sort the data according to K and not O
        # This is accomplished by zipping the two lists in to a single list of tuples
        # Then sort over the first element (default) which is K value
        # Then unzip to get back two sorted lists

        # The returned list may be empty in which case unzipping will fail
        if len(K):

            zipped = list(zip(K, O))
            zipped.sort()

            K, O = zip(*zipped)

        plt.plot(K, O, 'k')

    axes = plt.gca()
    axes.set_xlim([0, 5])
    axes.set_ylim([0, 4.5])

    plt.xlabel(r"$\beta_c = k_\perp \rho_c$")
    plt.ylabel(r"$\omega / \omega_c$", rotation=0)
    plt.title("Roots of Dispersion Relation")


def main(D_omega: ("Plot D(omega) with k_perp = 1 fixed", "flag", "o")):

    if D_omega:
        plot_D_omega()

    # Default option
    else:
        # plot_D_omega()
        plot_D()
        # plot_D_roots()
        # plot_c_array(c_Gamma_array, start=0.1, end=4, samples=100)
        # plot_hypergeom()

    # plt.legend()
    plt.grid(True)

    plt.show()


c_2F3 = c_functions.hyp2F3
c_2F3.restype = ctypes.c_double
c_2F3.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)


def plot_hypergeom():

    xs = numpy.linspace(0, 10, 100, endpoint=False)
    size = xs.size

    ys = numpy.empty((size,))
    cys = numpy.empty((size,))

    for i in range(size):

        x = xs[i]

        ys[i] = hyp2f3(1, 0.5, -3.5, 2.5, -0.5, x)
        cys[i] = c_2F3(1, 0.5, -3.5, 2.5, -0.5, x)

    plt.plot(xs, ys)
    plt.plot(xs, cys, 'k.')
    plt.title('2F3 as a function of x')


    plt.figure()
    b2s = numpy.linspace(1, 11, 100, endpoint=False)

    for i in range(size):

        b2 = b2s[i]

        ys[i] = hyp2f3(1, 0.5, -3.5, b2, -0.5, 1.5)
        cys[i] = c_2F3(1, 0.5, -3.5, b2, -0.5, 1.5)

    plt.plot(b2s, ys)
    plt.plot(b2s, cys, 'k.')
    plt.title('2F3 as a function of b2')


    plt.figure()
    # Needs to start at 1.1 because hyp1f2 is undefined at non-positive integers
    b3s = numpy.linspace(1.1, -9, 100, endpoint=False)

    for i in range(size):

        b3 = b3s[i]

        ys[i] = hyp2f3(1, 0.5, -3.5, 2.5, b3, 1.5)
        cys[i] = c_2F3(1, 0.5, -3.5, 2.5, b3, 1.5)

    plt.plot(b3s, ys)
    plt.plot(b3s, cys, 'k.')
    plt.title("2F3 as a function of b3")


if __name__ == '__main__':
    plac.call(main)
