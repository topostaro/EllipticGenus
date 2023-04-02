r"""
Module contains common functions used in intermediate calculations
"""

import math
from sage.all import (
    PolynomialRing,
    LaurentPolynomialRing,
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

    - ``x`` -- the variable of the todd genus for a line bundle.

    - ``m`` -- an integer.

    OUTPUT:

    the truncation of todd class `\frac{x}{1 - e^{-x}}` at degree ``m``.

    EXAMPLE:

        sage: from sage.WeakJacobiForm.elliptic_genus.utils import todd_cut
        sage: todd_cut(_Z, 5)

    """
    t = _L(_z / (1 - exp(-_z)), degree=m + 1).polynomial()
    return t(x)


def exp_cut(x, m):  # mでカットオフ
    e = _L(exp(_z), degree=m + 1).polynomial()
    return e(x)


def cutoff_for_coeff(m):  # xの式に対するcutoff
    return lambda f: sum(c * mono for c, mono in f if mono.total_degree() <= m)


def cutoff(f, m):  # x, y, qの式に対して, xの上の字数をcutoffする
    return f.map_coefficients(lambda g: g.map_coefficients(cutoff_for_coeff(m)))


def homogeneous_for_coeff(m):  # xの式に対するhomogeneous part
    return lambda f: sum(c * mono for c, mono in f if mono.total_degree() == m)


def homogeneous_part(f, m):  # x, y, qの式に対して, xの上の字数をhomogeneous part
    return f.map_coefficients(lambda g: g.map_coefficients(homogeneous_for_coeff(m)))


_m = SymmetricFunctions(QQ).m()
_e = SymmetricFunctions(QQ).e()
# partition [i_1, .., i_n]をmonomial symmetric polynomialの指数と捉えて、それを基本対称式に変換する関数
def chernnum_from_partition(dim: int, part):
    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    def monomial(degs):
        result = 1
        for deg in degs:
            if (deg > dim) or deg == 1:
                return 0
            else:
                result = result * c[deg]
        return result

    ls = list(_e(_m(part)))
    return sum(c * monomial(degs) for degs, c in ls)
