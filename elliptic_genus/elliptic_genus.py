from sage.all import (
    PolynomialRing,
    LaurentPolynomialRing,
    LazyLaurentSeriesRing,
    QQ,
    Partitions,
    prod,
)
from utils import *

# qのk次までに必要な係数の計算
def ell_factor_coeff_degreewise(dim: int, k: int) -> list:
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y_")
    y_ = R1.gen()
    R = LazyLaurentSeriesRing(R1, "q_")
    q_ = R.gen()

    S0 = PolynomialRing(QQ, "c", dim + 1)
    c = S0.gens()  # Chern根の変数

    S1 = LaurentPolynomialRing(S0, "y_")
    y = S1.gen()

    def factor(x, m):  # 無限積の因子
        numer = (1 - y_ * q_ * exp_cut(-x, m)) * (1 - y_ ^ (-1) * q_ * exp_cut(x, m))
        denom = (1 - q_ * exp_cut(-x, m)) * (1 - q_ * exp_cut(x, m))

        return numer / denom

    # 変数xに関するqがk次までの式(上の式のnはi)
    def elliptic_factor(x, n):
        term1 = prod(
            cutoff(factor(x, dim), dim)(q_ ^ i) for i in range(1, n + 1)
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
    s = LazyLaurentSeriesRing(S1, "q")
    q = S.gen()

    coeffs = [ell_factor_coeff_degreewise(dim, i) for i in range(k + 1)]

    return [sum(coeffs[i][j] * q**i for i in range(k + 1)) for j in range(dim + 1)]


def ell_coeff(dim: int, k: int) -> dict:
    R0 = PolynomialRing(QQ, "x", dim)
    x = R0.gens()  # Chern根の変数

    R1 = LaurentPolynomialRing(R0, "y")
    y = R1.gen()
    R = LazyLaurentSeriesRing(R1, "q")
    q = R.gen()

    efc = ell_factor_coeff(dim, k)

    return {
        part: (
            prod(efc[part[i]] for i in range(len(part))) * efc[0] ^ (dim - len(part))
        ).approximate_series(k + 1)
        for part in Partitions(dim)
    }


def elliptic_genus(dim: int, k: int):
    coeff = ell_coeff(dim, k)
    return sum(
        coeff[part] * chernnum_from_partition(dim, part) for part in Partitions(dim)
    )
