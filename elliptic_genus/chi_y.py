r"""
Computation of chi_y genus for complex manifolds
================================================

This module implements a computation of Eisenstein series.

EXAMPLES::


AUTHORS:

- KENTA KOBAYASHI (2023-03-30): initial version

REFERENCES:

"""

from sage.all import (
    PolynomialRing,
    LaurentPolynomialRing,
    LazyLaurentSeriesRing,
    QQ,
    Partitions,
    prod,
)
from utils import *

# chern rootごとのchi_yの積因子の係数
def _chi_y_factor_coeff(dim: int) -> list:
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y_")
    y_ = R1.gen()

    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    S1 = LaurentPolynomialRing(S0, "y")
    y = S1.gen()

    def chi_y_factor(x):
        term2 = 1 - y_ * exp_cut(-x, dim)
        term3 = todd_cut(x, dim)

        return (term2 * term3).map_coefficients(cutoff_for_coeff(dim))

    cy_fac_0 = chi_y_factor(x[0])
    result = [
        QQ(cy_fac_0.coefficients()[0].constant_coefficient())
        + QQ(cy_fac_0.coefficients()[1].constant_coefficient()) * y
    ] + [
        QQ(cy_fac_0.coefficients()[0].coefficient(x[0] ** i))
        + QQ(cy_fac_0.coefficients()[1].coefficient(x[0] ** i)) * y
        for i in range(1, dim + 1)
    ]
    return result


def _chi_y_coeff(dim: int) -> dict:
    cy_fac_coeff = _chi_y_factor_coeff(dim)
    return {
        part: prod(cy_fac_coeff[part[i]] for i in range(len(part)))
        * cy_fac_coeff[0] ** (dim - len(part))
        for part in Partitions(dim)
    }


def chi_y(dim: int):
    coeff = _chi_y_coeff(dim)
    return sum(
        coeff[part] * chernnum_from_partition(dim, part) for part in Partitions(dim)
    )


print(chi_y(3))
