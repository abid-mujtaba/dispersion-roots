"""
Test the ability to call C functions from Python.

The C functions are contained in a shared library: libfunctions.so

Source: https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
"""

import ctypes

_cfunctions = ctypes.CDLL('./libfunctions.so')          # Load the library

_cfunctions.J0.argtypes = (ctypes.c_double,)        # Explicitly declare the types of the arguments to the C function
_cfunctions.J0.restype = ctypes.c_double            # Explicity declare the type of the result of the function

_cfunctions.sum.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.sum.restype = ctypes.c_double

_cfunctions.triple.argypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int)
_cfunctions.triple.restype = None


def J0(x):

    ans = _cfunctions.J0(ctypes.c_double(x))        # The python object 'x' is cast to a ctypes.c_double type before passing it

    return float(ans)           # The result is explicitly cast back to python's 'float' type


def csum(xs):

    num = len(xs)               # C functions require the length of arrays to work
    array_type = ctypes.c_double * num      # Create a new type for an array of type c_double (with the length specified)

    # Note how we expand the list xs using * before passing it to array_type
    return float(_cfunctions.sum(array_type(*xs), ctypes.c_int(num)))


def triple(xs):

    num = len(xs)
    array_type = ctypes.c_double * num      # Create a new type for a double array of the specified length

    c_ys = array_type()     # Create a new C-type array that has the capacity to carry the result back
                            # Note: We didn't specify any list as the argument so an empty C-type array is created

    _cfunctions.triple(array_type(*xs), c_ys, num)      # Note how *xs needs to be cast to ctype but num can be sent as is and is automatically cast

    return [float(y) for y in c_ys]


def main():

    x = 5.0
    print("\nJ0({}) = {}".format(x, J0(x)))

    xs = [1.1, 2.2, 3.3, 4.4]
    print("Sum of {} = {}".format(str(xs), csum(xs)))

    ys = triple(xs)
    print("Tripled: {}".format(str(ys)))


if __name__ == '__main__':
    main()
