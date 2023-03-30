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
    r"""
    Return a function which is proposal to the Jacobi theta function $\theta_1$.

    $$
    \theta_1(\tau, z)
    = -\theta_{11}(\tau, z)
    = -i y^{1/2}q^{1/8}\sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2(n^2 + n)}
    $$
    $$
    t_1(y, q) = \sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2(n^2 + n)}
    $$
    where
    $q = e^{2 \pi i \tau}, y = e^{2 \pi i z}$

    INPUT:
    - `y` -- variable of the function which should be in `LaurentPolynomialRing(QQ, "y")`.

    OUTPUT:
    the function t_1 defined above.
    """

    # Return the pair of solutions of the equation "x^2 + x == n" or None.
    def sol(n):
        for i in range(n + 1):
            if i + i**2 == n:
                return (-i - 1, i)
        return None

    # Return the n-th coefficient.
    def coeff(n):
        s = sol(2 * n)
        if s == None:
            return 0
        else:
            s1, s2 = s
            return (-1) ** s1 * y**s1 + (-1) ** s2 * y**s2

    return _R(coeff, valuation=0)


def _function_t2(y):
    r"""
    Return a function which is proposal to the Jacobi theta function $\theta_2$.

    $$
    \theta_2(\tau, z)
    = \theta_{10}(\tau, z)
    = y^{1/2}q^{1/8}\sum_{n = -\infty}^\infty y^n q^{1/2(n^2 + n)}
    $$

    $$
    t_2(y, q) = \sum_{n = -\infty}^\infty y^n q^{1/2(n^2 + n)}
    $$

    INPUT:
    - `y` -- variable of the function which should be in `LaurentPolynomialRing(QQ, "y")`.

    OUTPUT:
    the function t_2 defined above.
    """

    # Return the pair of solutions of the equation "x^2 + x == n" or None.
    def sol(n):
        for i in range(n + 1):
            if i + i**2 == n:
                return (-i - 1, i)
        return None

    # Return the n-th coefficient.
    def coeff(n):
        s = sol(2 * n)
        if s == None:
            return 0
        else:
            s1, s2 = s
            return y**s1 + y**s2

    return _R(coeff, valuation=0)


def _function_t3(y):
    r"""
    Return a function which is proposal to the Jacobi theta function $\theta_3$.

    $$
    \theta_3(\tau, z)
    = \theta_{10}(\tau, z)
    = \sum_{n = -\infty}^\infty y^n q^{1/2n^2}
    $$

    $$
    t_3(y, q) = \sum_{n = -\infty}^\infty y^n q^{n^2}
    $$

    INPUT:
    - `y` -- variable of the function which should be in `LaurentPolynomialRing(QQ, "y")`.

    OUTPUT:
    the function t_3 defined above.
    """
    return _R(
        lambda n: 1
        if n == 0
        else y ** math.sqrt(n) + y ** (-math.sqrt(n))
        if math.floor(math.sqrt(n)) ** 2 == n
        else 0,
        valuation=0,
    )


def _function_t4(y):
    r"""
    Return a function which is proposal to the Jacobi theta function $\theta_3$.

    $$
    \theta_4(\tau, z)
    = \theta_{10}(\tau, z)
    = \sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2n^2}
    $$

    $$
    t_4(y, q) = \sum_{n = -\infty}^\infty (-1)^n y^n q^{n^2}
    $$

    INPUT:
    - `y` -- variable of the function which should be in `LaurentPolynomialRing(QQ, "y")`.

    OUTPUT:
    the function t_4 defined above.
    """
    return _R(
        lambda n: 1
        if n == 0
        else (-1) ** n * y ** math.sqrt(n) + (-1) ** n * y ** (-math.sqrt(n))
        if math.floor(math.sqrt(n)) ** 2 == n
        else 0,
        valuation=0,
    )


def _phi_tilde_0_1():
    r"""
    Return a weak Jacobi form of weight 0 and index 1.
    """
    f = 4 * (_function_t2(_y) / _function_t2(1)) ** 2 * _y

    g = (
        4 * (_function_t3(_y) / _function_t3(1)) ** 2
        + 4 * (_function_t4(_y) / _function_t4(1)) ** 2
    )
    return f + _R(lambda n: g.coefficient(2 * n), valuation=0)


def _phi_tilde_m2_1():
    r"""
    Return a weak Jacobi form of weight -2 and index 1.
    """
    h = _function_t3(1) * _function_t4(1)
    e3 = 1 / 2 * _function_t2(1) * _R(lambda n: h.coefficient(2 * n), valuation=0)
    return _function_t1(_y) ** 2 / e3**2 * _y


def _phi_tilde_0_3o2():
    r"""
    Return a weak Jacobi form of weight 0 and index 3/2 multiplied by $y^{1/2}$.
    """
    return _function_t1(_y**2) / _function_t1(_y) * _y


# nを2と3に分解するプログラム
def _decompose(n):
    r"""
    Return the list of all possible decompositions of `n` into sums of 2 and 3.

    INPUT:
    - `n` -- an integer expected to be positive for meaningful calculations.

    OUTPUT:
    The list of pairs of coefficient of possible decompositions of `n` into sums of 2 and 3.

    EXAMPLES::

        sage: _decompose(4)
        {(2, 0)}
        sage: _decompose(8)
        {(4, 0), (1, 2)}
    """
    if n <= 0:
        return []
    if n == 2:
        return [(1, 0)]
    if n == 3:
        return [(0, 1)]

    with2 = [(i + 1, j) for i, j in _decompose(n - 2)]
    with3 = [(i, j + 1) for i, j in _decompose(n - 3)]

    total = set(with2 + with3)
    return total


# 整数indexの基底を返す関数
def basis_integral(index: int) -> list:
    r"""
    Return a basis of the space of weak Jacobi forms of weight `0` and index `index`.

    INPUT:
    - `index` -- integral index of weak Jacobi forms expected to be positive.

    OUTPUT:
    the list consisting of series in a basis of weight 0 and index `index`.

    EXAMPLES::

        sage: basis(3)
        [(y^-3 + 30*y^-2 + 303*y^-1 + 1060 + 303*y + 30*y^2 + y^3) + (30*y^-4 + 408*y^-3 - 456*y^-2 - 12696*y^-1 + 25428 - 12696*y - 456*y^2 + 408*y^3 + 30*y^4)*q + (303*y^-5 - 456*y^-4 - 14085*y^-3 + 143280*y^-2 - 477738*y^-1 + 697392 - 477738*y + 143280*y^2 - 14085*y^3 - 456*y^4 + 303*y^5)*q^2 + (1060*y^-6 - 12696*y^-5 + 143280*y^-4 - 1070296*y^-3 + 4380348*y^-2 - 9943440*y^-1 + 13003488 - 9943440*y + 4380348*y^2 - 1070296*y^3 + 143280*y^4 - 12696*y^5 + 1060*y^6)*q^3 + (303*y^-7 + 25428*y^-6 - 477738*y^-5 + 4380348*y^-4 - 22852404*y^-3 + 71497068*y^-2 - 139461585*y^-1 + 173777160 - 139461585*y + 71497068*y^2 - 22852404*y^3 + 4380348*y^4 - 477738*y^5 + 25428*y^6 + 303*y^7)*q^4 + (30*y^-8 - 12696*y^-7 + 697392*y^-6 - 9943440*y^-5 + 71497068*y^-4 - 304558128*y^-3 + 825547488*y^-2 - 1481082024*y^-1 + 1795708620 - 1481082024*y + 825547488*y^2 - 304558128*y^3 + 71497068*y^4 - 9943440*y^5 + 697392*y^6 - 12696*y^7 + 30*y^8)*q^5 + (y^-9 - 456*y^-8 - 477738*y^-7 + 13003488*y^-6 - 139461585*y^-5 + 825547488*y^-4 - 3056436603*y^-3 + 7532417376*y^-2 - 12788771979*y^-1 + 15228360016 - 12788771979*y + 7532417376*y^2 - 3056436603*y^3 + 825547488*y^4 - 139461585*y^5 + 13003488*y^6 - 477738*y^7 - 456*y^8 + y^9)*q^6 + O(q^7), (y^-3 + 6*y^-2 - 33*y^-1 + 52 - 33*y + 6*y^2 + y^3) + (6*y^-4 + 120*y^-3 + 2040*y^-2 - 9336*y^-1 + 14340 - 9336*y + 2040*y^2 + 120*y^3 + 6*y^4)*q + (-33*y^-5 + 2040*y^-4 - 30501*y^-3 + 169776*y^-2 - 436410*y^-1 + 590256 - 436410*y + 169776*y^2 - 30501*y^3 + 2040*y^4 - 33*y^5)*q^2 + (52*y^-6 - 9336*y^-5 + 169776*y^-4 - 1238200*y^-3 + 4575660*y^-2 - 9643728*y^-1 + 12291552 - 9643728*y + 4575660*y^2 - 1238200*y^3 + 169776*y^4 - 9336*y^5 + 52*y^6)*q^3 + (-33*y^-7 + 14340*y^-6 - 436410*y^-5 + 4575660*y^-4 - 23929812*y^-3 + 72597180*y^-2 - 137835345*y^-1 + 170028840 - 137835345*y + 72597180*y^2 - 23929812*y^3 + 4575660*y^4 - 436410*y^5 + 14340*y^6 - 33*y^7)*q^4 + (6*y^-8 - 9336*y^-7 + 590256*y^-6 - 9643728*y^-5 + 72597180*y^-4 - 310041072*y^-3 + 830703072*y^-2 - 1473696456*y^-1 + 1779000156 - 1473696456*y + 830703072*y^2 - 310041072*y^3 + 72597180*y^4 - 9643728*y^5 + 590256*y^6 - 9336*y^7 + 6*y^8)*q^5 + (y^-9 + 2040*y^-8 - 436410*y^-7 + 12291552*y^-6 - 137835345*y^-5 + 830703072*y^-4 - 3080138715*y^-3 + 7553501280*y^-2 - 12759364635*y^-1 + 15162554320 - 12759364635*y + 7553501280*y^2 - 3080138715*y^3 + 830703072*y^4 - 137835345*y^5 + 12291552*y^6 - 436410*y^7 + 2040*y^8 + y^9)*q^6 + O(q^7), (y^-3 - 6*y^-2 + 15*y^-1 - 20 + 15*y - 6*y^2 + y^3) + (-6*y^-4 - 456*y^-3 + 2856*y^-2 - 7224*y^-1 + 9660 - 7224*y + 2856*y^2 - 456*y^3 - 6*y^4)*q + (15*y^-5 + 2856*y^-4 - 40005*y^-3 + 182160*y^-2 - 414666*y^-1 + 539280 - 414666*y + 182160*y^2 - 40005*y^3 + 2856*y^4 + 15*y^5)*q^2 + (-20*y^-6 - 7224*y^-5 + 182160*y^-4 - 1325176*y^-3 + 4671156*y^-2 - 9491280*y^-1 + 11940768 - 9491280*y + 4671156*y^2 - 1325176*y^3 + 182160*y^4 - 7224*y^5 - 20*y^6)*q^3 + (15*y^-7 + 9660*y^-6 - 414666*y^-5 + 4671156*y^-4 - 24474996*y^-3 + 73142916*y^-2 - 137017041*y^-1 + 168165912 - 137017041*y + 73142916*y^2 - 24474996*y^3 + 4671156*y^4 - 414666*y^5 + 9660*y^6 + 15*y^7)*q^4 + (-6*y^-8 - 7224*y^-7 + 539280*y^-6 - 9491280*y^-5 + 73142916*y^-4 - 312795504*y^-3 + 833272224*y^-2 - 1469993736*y^-1 + 1770666660 - 1469993736*y + 833272224*y^2 - 312795504*y^3 + 73142916*y^4 - 9491280*y^5 + 539280*y^6 - 7224*y^7 - 6*y^8)*q^5 + (y^-9 + 2856*y^-8 - 414666*y^-7 + 11940768*y^-6 - 137017041*y^-5 + 833272224*y^-4 - 3092014395*y^-3 + 7564027680*y^-2 - 12744642603*y^-1 + 15129690352 - 12744642603*y + 7564027680*y^2 - 3092014395*y^3 + 833272224*y^4 - 137017041*y^5 + 11940768*y^6 - 414666*y^7 + 2856*y^8 + y^9)*q^6 + O(q^7)]

    """
    result = []
    result.append(_phi_tilde_0_1() ** index)
    for i in range(2, index + 1):
        phi = _phi_tilde_0_1() ** (index - i) * _phi_tilde_m2_1() ** i
        for i, j in _decompose(i):
            e = eisenstein(4) ** i * eisenstein(6) ** j
            result.append(e * phi)

    return result


def basis_half_integral(double_index: int) -> list:
    r"""
    Return a basis of the space of weak Jacobi forms of weight `0` and index `double_index / 2` for even `double_index`,
    otherwise return a list of a basis multiplied by by $y^{1/2}$.

    INPUT:
    - `double_index` -- double of index of weak Jacobi forms expected to be greater than 3.

    OUTPUT:
    If `double_index` is even, the list consisting of series in a basis of weight `0` and index `double_index`.
    Otherwise, the list consisting of series in a basis of weight `0` and index `double_index` multiplied by $y^{1/2}$.
    """
    if double_index % 2 == 0:
        return basis_integral(int(double_index / 2))
    else:
        index = int((double_index - 3) / 2)
        return list(map(lambda f: f * _phi_tilde_0_3o2(), basis_integral(index)))
