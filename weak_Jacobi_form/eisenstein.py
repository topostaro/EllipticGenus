from sage.all import LaurentPolynomialRing, LazyLaurentSeriesRing, QQ, bernoulli, sigma

# Base rings for calculations
R0 = LaurentPolynomialRing(QQ, "y")
R = LazyLaurentSeriesRing(R0, "q")

# modified version only for negative integers
def _zeta(s: int):
    r"""
    Return the special values of zeta function only for negative integers
    """
    if s >= 0:
        raise TypeError("The argument must be negative.")
    return -bernoulli(-s + 1) / (-s + 1)


# Compute the coefficient of q^n
def _eisenstein_coefficient(k: int, n: int):
    r"""
    Return the n-th coefficient of Eisenstein series of index k.
    """
    if n == 0:
        return 1
    else:
        return 2 * sigma(n, k - 1) / _zeta(1 - k)


def eisenstein(k: int):
    r"""
    Return Eisenstein series of index k.

    EXAMPLES::

        sage: eisenstein(4)
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + O(q^7)
    """
    return R(lambda n: _eisenstein_coefficient(k, n), valuation=0)


print(eisenstein(4))
