r"""
Computation of elliptic genera of manifolds by using Chern numbers
================================================

This module implements a computation of elliptic genera, where the  coefficients are expressed by Chern numbers.

EXAMPLES:

    sage: from sage.EllipticGenus.elliptic_genus.elliptic_genus import elliptic_genus_chernnum
    sage: elliptic_genus_chernnum(4, 0)
    (-1/720*c1^4 + 1/180*c1^2*c2 + 1/240*c2^2 + 1/720*c1*c3 - 1/720*c4) + (1/180*c1^4 - 1/45*c1^2*c2 - 1/60*c2^2 + 7/90*c1*c3 + 31/180*c4)*y + (-1/120*c1^4 + 1/30*c1^2*c2 + 1/40*c2^2 - 19/120*c1*c3 + 79/120*c4)*y^2 + (1/180*c1^4 - 1/45*c1^2*c2 - 1/60*c2^2 + 7/90*c1*c3 + 31/180*c4)*y^3 + (-1/720*c1^4 + 1/180*c1^2*c2 + 1/240*c2^2 + 1/720*c1*c3 - 1/720*c4)*y^4 + O(q)
    sage: elliptic_genus_chernnum(5, 2)
    (-1/1440*c1^3*c2 + 1/480*c1*c2^2 + 1/1440*c1^2*c3 - 1/1440*c1*c4) + (1/480*c1^3*c2 - 1/160*c1*c2^2 - 1/480*c1^2*c3 + 7/160*c1*c4 + 1/24*c5)*y + (-1/720*c1^3*c2 + 1/240*c1*c2^2 + 1/720*c1^2*c3 - 31/720*c1*c4 + 11/24*c5)*y^2 + (-1/720*c1^3*c2 + 1/240*c1*c2^2 + 1/720*c1^2*c3 - 31/720*c1*c4 + 11/24*c5)*y^3 + (1/480*c1^3*c2 - 1/160*c1*c2^2 - 1/480*c1^2*c3 + 7/160*c1*c4 + 1/24*c5)*y^4 + (-1/1440*c1^3*c2 + 1/480*c1*c2^2 + 1/1440*c1^2*c3 - 1/1440*c1*c4)*y^5 + ((-1/24*c1^5 + 187/1440*c1^3*c2 - 7/480*c1*c2^2 - 247/1440*c1^2*c3 + 187/1440*c1*c4 - 1/24*c5)*y^-1 + (5/24*c1^5 - 235/288*c1^3*c2 + 55/96*c1*c2^2 + 151/288*c1^2*c3 - 295/288*c1*c4) + (-3/8*c1^5 + 267/160*c1^3*c2 - 261/160*c1*c2^2 - 87/160*c1^2*c3 + 207/160*c1*c4 - 9/4*c5)*y + (5/24*c1^5 - 283/288*c1^3*c2 + 103/96*c1*c2^2 + 55/288*c1^2*c3 - 115/288*c1*c4 + 55/24*c5)*y^2 + (5/24*c1^5 - 283/288*c1^3*c2 + 103/96*c1*c2^2 + 55/288*c1^2*c3 - 115/288*c1*c4 + 55/24*c5)*y^3 + (-3/8*c1^5 + 267/160*c1^3*c2 - 261/160*c1*c2^2 - 87/160*c1^2*c3 + 207/160*c1*c4 - 9/4*c5)*y^4 + (5/24*c1^5 - 235/288*c1^3*c2 + 55/96*c1*c2^2 + 151/288*c1^2*c3 - 295/288*c1*c4)*y^5 + (-1/24*c1^5 + 187/1440*c1^3*c2 - 7/480*c1*c2^2 - 247/1440*c1^2*c3 + 187/1440*c1*c4 - 1/24*c5)*y^6)*q + ((1/6*c1^5 - 83/240*c1^3*c2 + 3/80*c1*c2^2 + 1/80*c1^2*c3 + 127/240*c1*c4 - 11/24*c5)*y^-2 + (-23/12*c1^5 + 3049/480*c1^3*c2 - 609/160*c1*c2^2 - 363/160*c1^2*c3 + 243/160*c1*c4 + 9/4*c5)*y^-1 + (85/12*c1^5 - 2695/96*c1^3*c2 + 735/32*c1*c2^2 + 309/32*c1^2*c3 - 405/32*c1*c4) + (-137/12*c1^5 + 23791/480*c1^3*c2 - 7431/160*c1*c2^2 - 2477/160*c1^2*c3 + 9871/480*c1*c4 - 197/12*c5)*y + (73/12*c1^5 - 13199/480*c1^3*c2 + 4359/160*c1*c2^2 + 1293/160*c1^2*c3 - 1593/160*c1*c4 + 117/8*c5)*y^2 + (73/12*c1^5 - 13199/480*c1^3*c2 + 4359/160*c1*c2^2 + 1293/160*c1^2*c3 - 1593/160*c1*c4 + 117/8*c5)*y^3 + (-137/12*c1^5 + 23791/480*c1^3*c2 - 7431/160*c1*c2^2 - 2477/160*c1^2*c3 + 9871/480*c1*c4 - 197/12*c5)*y^4 + (85/12*c1^5 - 2695/96*c1^3*c2 + 735/32*c1*c2^2 + 309/32*c1^2*c3 - 405/32*c1*c4)*y^5 + (-23/12*c1^5 + 3049/480*c1^3*c2 - 609/160*c1*c2^2 - 363/160*c1^2*c3 + 243/160*c1*c4 + 9/4*c5)*y^6 + (1/6*c1^5 - 83/240*c1^3*c2 + 3/80*c1*c2^2 + 1/80*c1^2*c3 + 127/240*c1*c4 - 11/24*c5)*y^7)*q^2 + O(q^3)

AUTHORS:

- KENTA KOBAYASHI (2023-04-03): initial version

REFERENCES:

.. [EZ1985]  \Martin Eichler and Don Zagier. The theory of Jacobi forms, volume 55 of Progress in Mathematics. Birkh äuser Boston, Inc., Boston, MA, 1985.

.. [BL2000] \Lev A. Borisov and Anatoly Libgober. Elliptic genera of toric varieties and applications to mirror symmetry. Invent. Math., 140(2):453-485, 2000.
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
    LaurentPolynomialRing,
    LazyLaurentSeriesRing,
    QQ,
    Partitions,
    prod,
)
from sage.EllipticGenus.elliptic_genus.utils import *

# qのk次までに必要な係数の計算
def ell_factor_coeff_degreewise(dim: int, k: int) -> list:
    r"""

    Return the list which n-th component is the coefficients of `x^n q^k` of the following expression.

    ..MATH::

    \prod_{n=1}^\infty \frac{(1 - yq^ne^{-x})(1 - y^{-1}q^ne^x)}{(1 - q^ne^{-x})(1 - q^ne^x)}
    \cdot
    (1 - ye^{-x})
    \cdot
    \frac{x}{1 - e^{-x}}

    INPUT:

    - ``dim`` -- integer -- the dimension of the considering manifold

    - ``k`` -- integer

    OUTPUT:

    the list which n-th component is the coefficients of `x^n q^k` of the above expression.

    EXAMPLE:

        sage: from sage.EllipticGenus.elliptic_genus.elliptic_genus import ell_factor_coeff_degreewise
        sage: ell_factor_coeff_degreewise(5, 2)
        [-3*y^-1 + 9 - 9*y + 3*y^2,
         -9/2*y^-1 + 9/2 + 9/2*y - 9/2*y^2,
         -17/4*y^-1 + 27/4 - 27/4*y + 17/4*y^2,
         -3*y^-1 + 3 + 3*y - 3*y^2,
         -133/80*y^-1 + 159/80 - 159/80*y + 133/80*y^2,
         -3/4*y^-1 + 3/4 + 3/4*y - 3/4*y^2]

    """
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y_")
    y_ = R1.gen()
    R = LazyLaurentSeriesRing(R1, "q_")
    q_ = R.gen()

    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    S1 = LaurentPolynomialRing(S0, "y")
    y = S1.gen()

    def factor(x, m):  # 無限積の因子
        numer = (1 - y_ * q_ * exp_cut(-x, m)) * (1 - y_ ** (-1) * q_ * exp_cut(x, m))
        denom = (1 - q_ * exp_cut(-x, m)) * (1 - q_ * exp_cut(x, m))

        return numer / denom

    # 変数xに関するqがk次までの式(上の式のnはi)
    def elliptic_factor(x, n):
        term1 = prod(
            cutoff(factor(x, dim), dim)(q_**i) for i in range(1, n + 1)
        ).approximate_series(n + 1)
        term2 = 1 - y_ * exp_cut(-x, dim)
        term3 = todd_cut(x, dim)

        return cutoff(term1 * term2 * term3, dim)

    ell_fac = elliptic_factor(x[0], k + 1)[k]
    # y多項式としてみた時の係数の配列
    ell_fac_yc = list(ell_fac)
    val = ell_fac.valuation()

    result_element = 0
    for j in range(len(ell_fac_yc)):
        result_element += QQ(ell_fac_yc[j].constant_coefficient()) * y ** (j + val)

    result = [result_element]

    for i in range(1, dim + 1):
        result_element = 0
        for j in range(len(ell_fac_yc)):
            result_element += QQ(ell_fac_yc[j].coefficient(x[0] ** i)) * y ** (j + val)

        result += [result_element]

    return result


# k次までの和
def ell_factor_coeff(dim: int, k: int) -> list:
    r"""

    Return the list which n-th component is the coefficients of `x^n` of the following expression after taking the truncation of the terms of q variable up to degree ``k``.

    ..MATH::

    \prod_{n=1}^\infty \frac{(1 - yq^ne^{-x})(1 - y^{-1}q^ne^x)}{(1 - q^ne^{-x})(1 - q^ne^x)}
    \cdot
    (1 - ye^{-x})
    \cdot
    \frac{x}{1 - e^{-x}}

    INPUT:

    - ``dim`` -- integer -- the dimension of the considering manifold

    - ``k`` -- integer

    OUTPUT:

    the list which n-th component is the coefficients of `x^n` of the following expression after taking the truncation of the terms of q variable up to degree ``k``.

    EXAMPLE:

        sage: from sage.EllipticGenus.elliptic_genus.elliptic_genus import ell_factor_coeff
        sage: ell_factor_coeff(5, 2)
        [(1 - y) + (-y^-1 + 3 - 3*y + y^2)*q + (-3*y^-1 + 9 - 9*y + 3*y^2)*q^2,
         (1/2 + 1/2*y) + (-3/2*y^-1 + 3/2 + 3/2*y - 3/2*y^2)*q + (-9/2*y^-1 + 9/2 + 9/2*y - 9/2*y^2)*q^2,
         (1/12 - 1/12*y) + (-13/12*y^-1 + 5/4 - 5/4*y + 13/12*y^2)*q + (-17/4*y^-1 + 27/4 - 27/4*y + 17/4*y^2)*q^2,
         (-1/2*y^-1 + 1/2 + 1/2*y - 1/2*y^2)*q + (-3*y^-1 + 3 + 3*y - 3*y^2)*q^2,
         (-1/720 + 1/720*y) + (-119/720*y^-1 + 13/80 - 13/80*y + 119/720*y^2)*q + (-133/80*y^-1 + 159/80 - 159/80*y + 133/80*y^2)*q^2,
         (-1/24*y^-1 + 1/24 + 1/24*y - 1/24*y^2)*q + (-3/4*y^-1 + 3/4 + 3/4*y - 3/4*y^2)*q^2]

    """
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y_")
    y_ = R1.gen()
    R = LazyLaurentSeriesRing(R1, "q_")
    q_ = R.gen()

    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    S1 = LaurentPolynomialRing(S0, "y")
    y = S1.gen()
    S = LazyLaurentSeriesRing(S1, "q")
    q = S.gen()

    coeffs = [ell_factor_coeff_degreewise(dim, i) for i in range(k + 1)]

    return [sum(coeffs[i][j] * q**i for i in range(k + 1)) for j in range(dim + 1)]


def ell_coeff(dim: int, k: int) -> dict:
    r"""

    Return the dictionary which keys are partitions of ``dim`` and values are the coefficient ``x`` variable monomial with the multidegree of the key.

    INPUT:

    - ``dim`` -- integer -- the dimension of the considering manifold

    - ``k`` -- integer

    OUTPUT:

    the dictionary which keys are partitions of ``dim`` and values are the coefficient ``x`` variable monomial with the multidegree of the key.

    EXAMPLE:

        sage: from sage.EllipticGenus.elliptic_genus.elliptic_genus import ell_coeff
        sage: ell_coeff(5, 2)
        {[1, 1, 1, 1, 1]: 1/32 + 5/32*y + 5/16*y^2 + 5/16*y^3 + 5/32*y^4 + 1/32*y^5 + (-15/32*y^-1 - 45/32 - 15/32*y + 75/32*y^2 + 75/32*y^3 - 15/32*y^4 - 45/32*y^5 - 15/32*y^6)*q + (45/16*y^-2 + 45/32*y^-1 - 495/32 - 405/32*y + 765/32*y^2 + 765/32*y^3 - 405/32*y^4 - 495/32*y^5 + 45/32*y^6 + 45/16*y^7)*q^2 + O(q^3),
         [2, 1, 1, 1]: 1/96 + 1/96*y - 1/48*y^2 - 1/48*y^3 + 1/96*y^4 + 1/96*y^5 + (-23/96*y^-1 - 1/96 + 15/32*y - 7/32*y^2 - 7/32*y^3 + 15/32*y^4 - 1/96*y^5 - 23/96*y^6)*q + (83/48*y^-2 - 113/32*y^-1 - 77/32 + 779/96*y - 125/32*y^2 - 125/32*y^3 + 779/96*y^4 - 77/32*y^5 - 113/32*y^6 + 83/48*y^7)*q^2 + O(q^3),
         [2, 2, 1]: 1/288 - 1/96*y + 1/144*y^2 + 1/144*y^3 - 1/96*y^4 + 1/288*y^5 + (-31/288*y^-1 + 107/288 - 15/32*y + 59/288*y^2 + 59/288*y^3 - 15/32*y^4 + 107/288*y^5 - 31/288*y^6)*q + (17/16*y^-2 - 155/32*y^-1 + 313/32 - 349/32*y + 157/32*y^2 + 157/32*y^3 - 349/32*y^4 + 313/32*y^5 - 155/32*y^6 + 17/16*y^7)*q^2 + O(q^3),
         [3, 1, 1]: (-1/8*y^-1 + 1/8 + 3/8*y - 3/8*y^2 - 3/8*y^3 + 3/8*y^4 + 1/8*y^5 - 1/8*y^6)*q + (y^-2 - 15/4*y^-1 + 3/4 + 41/4*y - 33/4*y^2 - 33/4*y^3 + 41/4*y^4 + 3/4*y^5 - 15/4*y^6 + y^7)*q^2 + O(q^3),
         [3, 2]: (-1/24*y^-1 + 5/24 - 3/8*y + 5/24*y^2 + 5/24*y^3 - 3/8*y^4 + 5/24*y^5 - 1/24*y^6)*q + (2/3*y^-2 - 47/12*y^-1 + 115/12 - 143/12*y + 67/12*y^2 + 67/12*y^3 - 143/12*y^4 + 115/12*y^5 - 47/12*y^6 + 2/3*y^7)*q^2 + O(q^3),
         [4, 1]: -1/1440 + 1/480*y - 1/720*y^2 - 1/720*y^3 + 1/480*y^4 - 1/1440*y^5 + (-113/1440*y^-1 + 65/288 - 33/160*y + 17/288*y^2 + 17/288*y^3 - 33/160*y^4 + 65/288*y^5 - 113/1440*y^6)*q + (39/80*y^-2 - 517/160*y^-1 + 235/32 - 1203/160*y + 467/160*y^2 + 467/160*y^3 - 1203/160*y^4 + 235/32*y^5 - 517/160*y^6 + 39/80*y^7)*q^2 + O(q^3),
         [5]: (-1/24*y^-1 + 5/24 - 3/8*y + 5/24*y^2 + 5/24*y^3 - 3/8*y^4 + 5/24*y^5 - 1/24*y^6)*q + (1/6*y^-2 - 23/12*y^-1 + 85/12 - 137/12*y + 73/12*y^2 + 73/12*y^3 - 137/12*y^4 + 85/12*y^5 - 23/12*y^6 + 1/6*y^7)*q^2 + O(q^3)}

    """
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y")
    y = R1.gen()

    R = LazyLaurentSeriesRing(R1, "q")
    q = R.gen()

    efc = ell_factor_coeff(dim, k)

    return {
        part: (
            prod(efc[part[i]] for i in range(len(part))) * efc[0] ** (dim - len(part))
        ).approximate_series(k + 1)
        for part in Partitions(dim)
    }


def elliptic_genus_chernnum(dim: int, k: int):
    r"""

    Return the elliptic genus of the manifold of dimension ``dim`` with the terms of q variable up to degree ``k``,
    where the  coefficients are expressed by Chern numbers.

    INPUT:

    - ``dim`` -- integer -- the dimension of the considering manifold

    - ``k`` -- integer

    OUTPUT:

    the elliptic genus of the manifold of dimension ``dim`` with the terms of q variable up to degree ``k``,
    where the  coefficients are expressed by Chern numbers.

    EXAMPLE:

        sage: from sage.EllipticGenus.elliptic_genus.elliptic_genus import elliptic_genus_chernnum
        sage: elliptic_genus_chernnum(5, 2)
        (-1/1440*c1^3*c2 + 1/480*c1*c2^2 + 1/1440*c1^2*c3 - 1/1440*c1*c4) + (1/480*c1^3*c2 - 1/160*c1*c2^2 - 1/480*c1^2*c3 + 7/160*c1*c4 + 1/24*c5)*y + (-1/720*c1^3*c2 + 1/240*c1*c2^2 + 1/720*c1^2*c3 - 31/720*c1*c4 + 11/24*c5)*y^2 + (-1/720*c1^3*c2 + 1/240*c1*c2^2 + 1/720*c1^2*c3 - 31/720*c1*c4 + 11/24*c5)*y^3 + (1/480*c1^3*c2 - 1/160*c1*c2^2 - 1/480*c1^2*c3 + 7/160*c1*c4 + 1/24*c5)*y^4 + (-1/1440*c1^3*c2 + 1/480*c1*c2^2 + 1/1440*c1^2*c3 - 1/1440*c1*c4)*y^5 + ((-1/24*c1^5 + 187/1440*c1^3*c2 - 7/480*c1*c2^2 - 247/1440*c1^2*c3 + 187/1440*c1*c4 - 1/24*c5)*y^-1 + (5/24*c1^5 - 235/288*c1^3*c2 + 55/96*c1*c2^2 + 151/288*c1^2*c3 - 295/288*c1*c4) + (-3/8*c1^5 + 267/160*c1^3*c2 - 261/160*c1*c2^2 - 87/160*c1^2*c3 + 207/160*c1*c4 - 9/4*c5)*y + (5/24*c1^5 - 283/288*c1^3*c2 + 103/96*c1*c2^2 + 55/288*c1^2*c3 - 115/288*c1*c4 + 55/24*c5)*y^2 + (5/24*c1^5 - 283/288*c1^3*c2 + 103/96*c1*c2^2 + 55/288*c1^2*c3 - 115/288*c1*c4 + 55/24*c5)*y^3 + (-3/8*c1^5 + 267/160*c1^3*c2 - 261/160*c1*c2^2 - 87/160*c1^2*c3 + 207/160*c1*c4 - 9/4*c5)*y^4 + (5/24*c1^5 - 235/288*c1^3*c2 + 55/96*c1*c2^2 + 151/288*c1^2*c3 - 295/288*c1*c4)*y^5 + (-1/24*c1^5 + 187/1440*c1^3*c2 - 7/480*c1*c2^2 - 247/1440*c1^2*c3 + 187/1440*c1*c4 - 1/24*c5)*y^6)*q + ((1/6*c1^5 - 83/240*c1^3*c2 + 3/80*c1*c2^2 + 1/80*c1^2*c3 + 127/240*c1*c4 - 11/24*c5)*y^-2 + (-23/12*c1^5 + 3049/480*c1^3*c2 - 609/160*c1*c2^2 - 363/160*c1^2*c3 + 243/160*c1*c4 + 9/4*c5)*y^-1 + (85/12*c1^5 - 2695/96*c1^3*c2 + 735/32*c1*c2^2 + 309/32*c1^2*c3 - 405/32*c1*c4) + (-137/12*c1^5 + 23791/480*c1^3*c2 - 7431/160*c1*c2^2 - 2477/160*c1^2*c3 + 9871/480*c1*c4 - 197/12*c5)*y + (73/12*c1^5 - 13199/480*c1^3*c2 + 4359/160*c1*c2^2 + 1293/160*c1^2*c3 - 1593/160*c1*c4 + 117/8*c5)*y^2 + (73/12*c1^5 - 13199/480*c1^3*c2 + 4359/160*c1*c2^2 + 1293/160*c1^2*c3 - 1593/160*c1*c4 + 117/8*c5)*y^3 + (-137/12*c1^5 + 23791/480*c1^3*c2 - 7431/160*c1*c2^2 - 2477/160*c1^2*c3 + 9871/480*c1*c4 - 197/12*c5)*y^4 + (85/12*c1^5 - 2695/96*c1^3*c2 + 735/32*c1*c2^2 + 309/32*c1^2*c3 - 405/32*c1*c4)*y^5 + (-23/12*c1^5 + 3049/480*c1^3*c2 - 609/160*c1*c2^2 - 363/160*c1^2*c3 + 243/160*c1*c4 + 9/4*c5)*y^6 + (1/6*c1^5 - 83/240*c1^3*c2 + 3/80*c1*c2^2 + 1/80*c1^2*c3 + 127/240*c1*c4 - 11/24*c5)*y^7)*q^2 + O(q^3)

    """

    coeff = ell_coeff(dim, k)
    return sum(
        coeff[part] * chernnum_from_partition(dim, part) for part in Partitions(dim)
    )
