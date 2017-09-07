from mpmath import hyp1f2, hyp2f3, gamma, csc
from math import sqrt, pi
import matplotlib.pyplot as plt
import numpy


def main():

    k = 1.0 
    rho_h = 1
    Kh = 4

    a1 = 1
    a2 = 0.5
    b1 = 0.5 - Kh

    two_lambda_j = 2 * (Kh - 1.5) * k**2 * rho_h**2
    
    # Define the 2F3 as a function of omega only the rest being constants we have defined.
    def F23(w):
        return hyp2f3(a1, a2, b1, 1 + w, 1 - w, two_lambda_j)

    # Plot values of F23(w)
    def plot():

        omega = numpy.linspace(2.1, 7.9, 10000)
        F = numpy.empty_like(omega)

        # Apply the function F23 to the values in the array 'omega'
        for i in range(omega.size):
            F[i] = F23( omega[i] )

        plt.plot(omega, F, '.') 
        #plt.ylim([-10,10])
        plt.grid()
        plt.show()


    # Print boundary and middle values of 2F3 in the bands
    def values():

        for i in range(1,8):
            w1 = i + 1e-6
            w2 = i + 0.5
            w3 = i + 1 - 1e-6

            for w in (w1, w2, w3):
                print("2F3({}) = {}".format(w, F23(w)))

            print()


    values()




def D(k, w):
    """
    Henning
    """

    two_lambda_j_prime = 2 * (4 - 1.5) * pow(k, 2)
    h2f3 = hyp2f3(1, 0.5, 0.5 - 4, 1 + w, 1 - w, two_lambda_j_prime)

    print(two_lambda_j_prime)
    print(int(h2f3))

    result = 1
    result -= h2f3

    print(result)

    third = sqrt(pi) * w
    third *= gamma(4 + 1) * gamma(0.5 - 4)
    third /= csc(pi * w)
    third /= gamma(4 + 1.5 + w) * gamma(4 + 1.5 - w)
    third *= pow(two_lambda_j_prime, 4 + 0.5)
    third *= hyp1f2(4 + 1, 4 + 1.5 + w, 4 + 1.5 - w, two_lambda_j_prime)

    print(int(third))

    result += third

    print(result)

    lambda_kappa_j_p2 = (4 - 1.5) / ((pow(5.099, 2) - 1) * (4 - 0.5))
    print(lambda_kappa_j_p2 * pow(k,2))

    result /= pow(k, 2) * lambda_kappa_j_p2

    return result


if __name__ == '__main__':
    main()
