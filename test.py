from mpmath import hyp1f2, hyp2f3, gamma, csc
from math import sqrt, pi


def main():

    k = 6.6
    w = 0.5

    print("D({:f}, {:f}) = {:f}".format(k, w, float(D(k, w))))


def D(k, w):

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
