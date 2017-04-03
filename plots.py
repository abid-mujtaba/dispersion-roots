"""
Create plots of the Dispersion Relation and its roots.

The python functions for calculating these values are implemented in functions.py
"""

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy
import plac

from functions import D_roots, D_mesh


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

    for i in range(1, 8):

        print("Calculating roots between omega = {} and {}".format(i, i + 1))

        # The function is undefined at integer and half-integer values so we boost all values by 1e-3 to avoid these
        # The third argument is the number of points to be created in the interval.
        # We are ignoring the end-point to keep it within the region
        slices = numpy.linspace(i + 1e-3, i + 1, 200, endpoint=False)

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
    axes.set_xlim([0, 100])
    axes.set_ylim([1, 8])

    plt.xlabel(r"$\beta_c = k_\perp \rho_c$")
    plt.ylabel(r"$\omega / \omega_c$", rotation=0)
    plt.title("Roots of Dispersion Relation")


def main(D_omega: ("Plot D(omega) with k_perp = 1 fixed", "flag", "o")):

    if D_omega:
        plot_D_omega()

    # Default option
    else:
        # plot_D()
        plot_D_roots()

    # plt.legend()
    plt.grid(True)

    plt.show()


if __name__ == '__main__':
    plac.call(main)
