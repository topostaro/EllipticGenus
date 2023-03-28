from sage.all import *

# modified version only for negative integers
def zeta(s: int):
    return -bernoulli(-s + 1) / (-s + 1)


# Compute the coefficient of q^n
def eisenstein_coefficient(k: int, n: int):
    if n == 0:
        return 1
    else:
        return 2 * sigma(n, k - 1) / zeta(1 - k)


def eisenstein(k: int):
    return R(lambda n: eisenstein_coefficient(k, n), valuation=0)
