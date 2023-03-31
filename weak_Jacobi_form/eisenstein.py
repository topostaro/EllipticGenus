from sage.all import LaurentPolynomialRing, LazyLaurentSeriesRing, QQ, bernoulli, sigma

# Base rings for calculations
R0 = LaurentPolynomialRing(QQ, "y")
R = LazyLaurentSeriesRing(R0, "q")

# modified version only for negative integers
def _zeta(s: int):
    r"""

    Return the special values of zeta function only for negative integers

    INPUT:

    - ``s`` -- a negative integer

    OUTPUT:

    The value of zeta function

    """
    if s >= 0:
        raise TypeError("The argument must be negative.")
    return -bernoulli(-s + 1) / (-s + 1)


# Compute the coefficient of q^n
def _eisenstein_coefficient(k: int, n: int):
    r"""

    Return the ``n``-th coefficient of Eisenstein series of index ``k``.

    INPUT:

    - ``k`` -- an integer representing the index of considering Eisenstein series.
    - ``n`` -- an integer representing the degree.

    OUTPUT:

    The ``n``-th coefficient of Eisenstein series of index ``k``.

    """
    if n == 0:
        return 1
    else:
        return 2 * sigma(n, k - 1) / _zeta(1 - k)


def eisenstein(k: int):
    r"""
    Return Eisenstein series of index ``k``.

    INPUT:

    - ``k`` -- an integer representing the index of considering Eisenstein series.

    OUTPUT:

    The Eisenstein series of index ``k``.

    EXAMPLES::

        sage: from sage.WeakJacobiForm.weak_Jacobi_form.eisenstein import eisenstein
        sage: eisenstein(4)
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + O(q^7)
        sage: eisenstein(6)
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 + O(q^7)
        sage: eisenstein(8)
        1 + 480*q + 61920*q^2 + 1050240*q^3 + 7926240*q^4 + 37500480*q^5 + 135480960*q^6 + O(q^7)

    """
    return R(lambda n: _eisenstein_coefficient(k, n), valuation=0)
