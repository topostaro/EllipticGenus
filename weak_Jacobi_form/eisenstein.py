r"""
Computation of Eisenstein series
================================================

This module implements a computation of Eisenstein series.

EXAMPLES::

    sage: from weak_Jacobi_form.eisenstein import eisenstein
    sage: lazy_eisenstein_series_qexp(4)
    1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + O(q^7)
    sage: lazy_eisenstein_series_qexp(6)
    1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 + O(q^7)
    sage: lazy_eisenstein_series_qexp(8, 5)
    1 + 480*q + 61920*q^2 + 1050240*q^3 + 7926240*q^4 + O(q^5)


AUTHORS:

- KENTA KOBAYASHI (2023-03-30): initial version

REFERENCES:

.. [EZ1985]  \Martin Eichler and Don Zagier. The theory of Jacobi forms, volume 55 of Progress in Mathematics. Birkh äuser Boston, Inc., Boston, MA, 1985.

"""

# ****************************************************************************
#       Copyright (C) 2023 KENTA KOBAYASHI <kenta.topos@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import LazyLaurentSeriesRing, QQ, bernoulli, sigma


# modified version only for negative integers
def _zeta(s: int):
    r"""

    Return the special values of zeta function only for negative integers

    INPUT:

    - ``s`` -- negative integer

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

    - ``k`` -- integer -- the index of considering Eisenstein series.
    - ``n`` -- integer -- the degree.

    OUTPUT:

    The ``n``-th coefficient of Eisenstein series of index ``k``.

    """
    if n == 0:
        return 1
    else:
        return 2 * sigma(n, k - 1) / _zeta(1 - k)


def lazy_eisenstein_series_qexp(k: int, prec=10, K=QQ, var="q"):
    r"""
    Return Eisenstein series of index ``k`` expressed by lazy series

    INPUT:

    - ``k`` -- integer -- the index of considering Eisenstein series.

    - ``prec`` -- integer -- (default: 10) the precision

    - ``K`` -- (default: `\QQ`) the base ring

    - ``var`` -- (default: ``'q'``) the variable name to use for q-expansion

    OUTPUT:

    The Eisenstein series of index ``k``  expressed by lazy series.

    EXAMPLES::

        sage: from weak_Jacobi_form.eisenstein import eisenstein
        sage: lazy_eisenstein_series_qexp(4)
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + 60480*q^6 + O(q^7)
        sage: lazy_eisenstein_series_qexp(6)
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 - 4058208*q^6 + O(q^7)
        sage: lazy_eisenstein_series_qexp(8, 5)
        1 + 480*q + 61920*q^2 + 1050240*q^3 + 7926240*q^4 + O(q^5)

    """

    # Base rings for calculations
    R = LazyLaurentSeriesRing(K, var)

    return R(lambda n: _eisenstein_coefficient(k, n), valuation=0).approximate_series(
        prec
    )
