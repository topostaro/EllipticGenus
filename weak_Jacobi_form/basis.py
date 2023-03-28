# from sage.rings.rational_field import QQ
# from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
# from sage.rings.lazy_series_ring import LazyLaurentSeriesRing

import math
from sage.all import LaurentPolynomialRing, LazyLaurentSeriesRing, QQ
from eisenstein import eisenstein

_R0 = LaurentPolynomialRing(QQ, "y")
_y = _R0.gen()

_R = LazyLaurentSeriesRing(_R0, "q")
_q = _R.gen()


def _function_t1(y):
    def sol(n):
        for i in range(n + 1):
            if i + i**2 == n:
                return (-i - 1, i)
        return None

    def coeff(n):
        s = sol(2 * n)
        if s == None:
            return 0
        else:
            s1, s2 = s
            return (-1) ** s1 * y**s1 + (-1) ** s2 * y**s2

    return _R(coeff, valuation=0)


def _function_t2(y):
    def sol(n):
        for i in range(n + 1):
            if i + i**2 == n:
                return (-i - 1, i)
        return None

    def coeff(n):
        s = sol(2 * n)
        if s == None:
            return 0
        else:
            s1, s2 = s
            return y**s1 + y**s2

    return _R(coeff, valuation=0)


def _function_t3(y):
    return _R(
        lambda n: 1
        if n == 0
        else y ** math.sqrt(n) + y ** (-math.sqrt(n))
        if math.floor(math.sqrt(n)) ** 2 == n
        else 0,
        valuation=0,
    )


def _function_t4(y):
    return _R(
        lambda n: 1
        if n == 0
        else (-1) ** n * y ** math.sqrt(n) + (-1) ** n * y ** (-math.sqrt(n))
        if math.floor(math.sqrt(n)) ** 2 == n
        else 0,
        valuation=0,
    )


f = 4 * (_function_t2(_y) / _function_t2(1)) ** 2 * _y
g = (
    4 * (_function_t3(_y) / _function_t3(1)) ** 2
    + 4 * (_function_t4(_y) / _function_t4(1)) ** 2
)
_phi_tilde_0_1 = f + _R(lambda n: g.coefficient(2 * n), valuation=0)

h = _function_t3(1) * _function_t4(1)
e3 = 1 / 2 * _function_t2(1) * _R(lambda n: h.coefficient(2 * n), valuation=0)
_phi_tilde_m2_1 = _function_t1(_y) ** 2 / e3**2 * _y


# nを2と3に分解するプログラム
def decompose(n):
    if n <= 0:
        return []
    if n == 2:
        return [(1, 0)]
    if n == 3:
        return [(0, 1)]

    with2 = [(i + 1, j) for i, j in decompose(n - 2)]
    with3 = [(i, j + 1) for i, j in decompose(n - 3)]

    total = set(with2 + with3)
    return total


# 整数indexの基底を返す関数
def basis(index: int) -> list:
    result = []
    result.append(_phi_tilde_0_1**index)
    for i in range(2, index + 1):
        phi = _phi_tilde_0_1 ** (index - i) * _phi_tilde_m2_1**i
        for i, j in decompose(i):
            e = eisenstein(4) ** i * eisenstein(6) ** j
            result.append(e * phi)

    return result


print(basis(3))
