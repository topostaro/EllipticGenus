r"""
Module contains common functions used in intermediate calculations
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

from sage.all import (
    PolynomialRing,
    LazyLaurentSeriesRing,
    QQ,
    SymmetricFunctions,
    exp,
)

_L = LazyLaurentSeriesRing(QQ, "z")
_z = _L.gen()


def todd_cut(x, m):  # mでカットオフ
    r"""
    Return the truncation of `\frac{x}{1 - e^{-x}}` at degree ``m``.

    INPUT:

    - ``x`` -- the variable of the Todd genus for a line bundle.

    - ``m`` -- integer.

    OUTPUT:

    the truncation of todd class `\frac{x}{1 - e^{-x}}` at degree ``m``.

    EXAMPLES::

        sage: from elliptic_genus.utils import todd_cut
        sage: R.<z> = LazyLaurentSeriesRing(QQ)
        sage: todd_cut(z, 5)
        1 + 1/2*z + 1/12*z^2 - 1/720*z^4

    """
    t = _L(_z / (1 - exp(-_z)), degree=m + 1).polynomial()
    return t(x)


def exp_cut(x, m):  # mでカットオフ
    r"""

    Return the truncation of `e^x` at degree ``m``.

    INPUT:

    - ``x`` -- the variable of the exponential function.

    - ``m`` -- integer.

    OUTPUT:

    the truncation of todd class `e^x` at degree ``m``.

    EXAMPLES::

        sage: from elliptic_genus.utils import exp_cut
        sage: R.<z> = LazyLaurentSeriesRing(QQ)
        sage: exp_cut(z, 5)
        1 + z + 1/2*z^2 + 1/6*z^3 + 1/24*z^4 + 1/120*z^5

    """
    e = _L(exp(_z), degree=m + 1).polynomial()
    return e(x)


def cutoff_for_coeff(m):  # xの式に対するcutoff
    r"""

    Return a function which returns the truncation of its argument
    at degree ``m``.

    INPUT:

    - ``m`` -- integer.

    OUTPUT:

    a function which returns the truncation of the argument function at degree ``m``.

    EXAMPLES::

        sage: from elliptic_genus.utils import cutoff_for_coeff
        sage: R0.<x0, x1> = PolynomialRing(QQ)
        sage: cutoff_for_coeff(2)(x0 + x0 * x1 - 2 * x1^4)
        x0*x1 + x0


    """
    return lambda f: sum(c * mono for c, mono in f if mono.total_degree() <= m)


def cutoff(f, m):  # x, y, qの式に対して, xの上の字数をcutoffする
    r"""

    Return the truncation of the argument truncating only the coefficients
    of the coefficients at degree``m``.

    INPUT:

    - ``f`` -- a polynomial with coefficients of polynomials with coefficients of polynomials.

    - ``m`` -- integer.

    OUTPUT:

    the truncation of the argument truncating only the coefficients
    of the coefficients.

    EXAMPLES::

        sage: from elliptic_genus.utils import cutoff
        sage: R0.<x0, x1> = PolynomialRing(QQ)
        sage: R1.<y> = LaurentPolynomialRing(R0)
        sage: R.<z>  = LazyLaurentSeriesRing(R1)
        sage: f = (1 + x0 + x1)^5 * y^6 * z^7
        sage: cutoff(f, 4)
        ((5*x0^4 + 20*x0^3*x1 + 30*x0^2*x1^2 + 20*x0*x1^3 + 5*x1^4 + 10*x0^3 + 30*x0^2*x1 + 30*x0*x1^2 + 10*x1^3 + 10*x0^2 + 20*x0*x1 + 10*x1^2 + 5*x0 + 5*x1 + 1)*y^6)*z^7

    """
    return f.map_coefficients(lambda g: g.map_coefficients(cutoff_for_coeff(m)))


def homogeneous_for_coeff(m):  # xの式に対するhomogeneous part
    r"""

    Return a function which returns the homogeneous part of its argument
    at degree ``m``.

    INPUT:

    - ``m`` -- integer.

    OUTPUT:

    a function which returns the homogeneous part of the argument function
    at degree ``m``.

    EXAMPLES::

        sage: from elliptic_genus.utils import homogeneous_for_coeff
        sage: R0.<x0, x1> = PolynomialRing(QQ)
        sage: homogeneous_for_coeff(2)(x0 + x0 * x1 - 2 * x1^4)
        x0*x1

    """
    return lambda f: sum(c * mono for c, mono in f if mono.total_degree() == m)


def homogeneous_part(f, m):  # x, y, qの式に対して, xの上の字数をhomogeneous part
    r"""

    Return the homogeneous part of the coefficients of the coefficients
    at degree ``m``.

    INPUT:

    - ``f`` -- a polynomial with coefficients of polynomials with coefficients of polynomials.

    - ``m`` -- integer.

    OUTPUT:

    the homogeneous part of the coefficients of the coefficients at degree ``m``.

    EXAMPLES::

        sage: from elliptic_genus.utils import homogeneous_part
        sage: R0.<x0, x1> = PolynomialRing(QQ)
        sage: R1.<y> = LaurentPolynomialRing(R0)
        sage: R.<z>  = LazyLaurentSeriesRing(R1)
        sage: f = (1 + x0 + x1)^5 * y^6 * z^7
        sage: homogeneous_part(f, 4)
        ((5*x0^4 + 20*x0^3*x1 + 30*x0^2*x1^2 + 20*x0*x1^3 + 5*x1^4)*y^6)*z^7

    """
    return f.map_coefficients(lambda g: g.map_coefficients(homogeneous_for_coeff(m)))


# partition [i_1, .., i_n]をmonomial symmetric polynomialの指数と捉えて、それを基本対称式に変換する関数
def chernnum_from_partition(dim: int, part):
    r"""

    Return the combination of elementary symmetric functions equals to
    the monomial expressed by the argument partition.

    This function is used to convert Chern roots' monomials into
    polynomials of Chern classes. Here, the monomial of degree ``dim``
    is represented as a partition ``part`` of ``dim``. Since Chern classes
    are equal to the elementary symmetric functions of Chern roots,
    we rewrite the monomial represented by ``part`` into the expression using
    the elementary symmetric functions and then replace
    the i-th elementary symmetric functions  with ``c_i``.

    INPUT:

    - ``dim`` -- integer -- the dimension of the considering manifold.

    - ``part`` -- partition of ``dim`` -- the multidegree of a monomial.

    OUTPUT:

    the combination of elementary symmetric functions equals to the monomial with multidegree ``part``.

    EXAMPLES::

        sage: from elliptic_genus.utils import chernnum_from_partition
        sage: chernnum_from_partition(5, [5])
        c1^5 - 5*c1^3*c2 + 5*c1*c2^2 + 5*c1^2*c3 - 5*c2*c3 - 5*c1*c4 + 5*c5
        sage: chernnum_from_partition(3, [2,1])
        c1*c2 - 3*c3

    """

    m = SymmetricFunctions(QQ).m()
    e = SymmetricFunctions(QQ).e()

    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    def monomial(degs):
        result = 1
        for deg in degs:
            if deg > dim:
                return 0
            else:
                result = result * c[deg]
        return result

    ls = list(e(m(part)))
    return sum(c * monomial(degs) for degs, c in ls)
