r"""
Computation of a basis of weak Jacobi forms of weight `0`
================================================

This module implements a computation of a basis of weak Jacobi forms of weight `0`
and indices that are either integral or half-integral.

EXAMPLES::

    sage: from sage.WeakJacobiForm.weak_Jacobi_form.basis import basis_integral, basis_half_integral
    sage: basis_integral(2)
    [(y^-2 + 20*y^-1 + 102 + 20*y + y^2) + (20*y^-3 + 72*y^-2 - 1044*y^-1 + 1904 - 1044*y + 72*y^2 + 20*y^3)*q + (102*y^-4 - 1044*y^-3 + 7392*y^-2 - 23532*y^-1 + 34164 - 23532*y + 7392*y^2 - 1044*y^3 + 102*y^4)*q^2 + (20*y^-5 + 1904*y^-4 - 23532*y^-3 + 115552*y^-2 - 283688*y^-1 + 379488 - 283688*y + 115552*y^2 - 23532*y^3 + 1904*y^4 + 20*y^5)*q^3 + (y^-6 - 1044*y^-5 + 34164*y^-4 - 283688*y^-3 + 1107711*y^-2 - 2402244*y^-1 + 3090200 - 2402244*y + 1107711*y^2 - 283688*y^3 + 34164*y^4 - 1044*y^5 + y^6)*q^4 + (72*y^-6 - 23532*y^-5 + 379488*y^-4 - 2402244*y^-3 + 8066472*y^-2 - 16135248*y^-1 + 20229984 - 16135248*y + 8066472*y^2 - 2402244*y^3 + 379488*y^4 - 23532*y^5 + 72*y^6)*q^5 + (20*y^-7 + 7392*y^-6 - 283688*y^-5 + 3090200*y^-4 - 16135248*y^-3 + 48552352*y^-2 - 91617180*y^-1 + 112772304 - 91617180*y + 48552352*y^2 - 16135248*y^3 + 3090200*y^4 - 283688*y^5 + 7392*y^6 + 20*y^7)*q^6 + O(q^7),
    (y^-2 - 4*y^-1 + 6 - 4*y + y^2) + (-4*y^-3 + 264*y^-2 - 1020*y^-1 + 1520 - 1020*y + 264*y^2 - 4*y^3)*q + (6*y^-4 - 1020*y^-3 + 8160*y^-2 - 23556*y^-1 + 32820 - 23556*y + 8160*y^2 - 1020*y^3 + 6*y^4)*q^2 + (-4*y^-5 + 1520*y^-4 - 23556*y^-3 + 117856*y^-2 - 283640*y^-1 + 375648 - 283640*y + 117856*y^2 - 23556*y^3 + 1520*y^4 - 4*y^5)*q^3 + (y^-6 - 1020*y^-5 + 32820*y^-4 - 283640*y^-3 + 1113855*y^-2 - 2402316*y^-1 + 3080600 - 2402316*y + 1113855*y^2 - 283640*y^3 + 32820*y^4 - 1020*y^5 + y^6)*q^4 + (264*y^-6 - 23556*y^-5 + 375648*y^-4 - 2402316*y^-3 + 8081256*y^-2 - 16135152*y^-1 + 20207712 - 16135152*y + 8081256*y^2 - 2402316*y^3 + 375648*y^4 - 23556*y^5 + 264*y^6)*q^5 + (-4*y^-7 + 8160*y^-6 - 283640*y^-5 + 3080600*y^-4 - 16135152*y^-3 + 48585376*y^-2 - 91617300*y^-1 + 112723920 - 91617300*y + 48585376*y^2 - 16135152*y^3 + 3080600*y^4 - 283640*y^5 + 8160*y^6 - 4*y^7)*q^6 + O(q^7)]
    sage: basis_half_integral(5)
    [(y^-1 + 11 + 11*y + y^2) + (-y^-3 - 54*y^-1 + 55 + 55*y - 54*y^2 - y^4)*q + (-11*y^-4 + 54*y^-3 - 394*y^-1 + 351 + 351*y - 394*y^2 + 54*y^4 - 11*y^5)*q^2 + (-11*y^-5 - 55*y^-4 + 394*y^-3 - 1889*y^-1 + 1561 + 1561*y - 1889*y^2 + 394*y^4 - 55*y^5 - 11*y^6)*q^3 + (-y^-6 - 55*y^-5 - 351*y^-4 + 1889*y^-3 - 7398*y^-1 + 5916 + 5916*y - 7398*y^2 + 1889*y^4 - 351*y^5 - 55*y^6 - y^7)*q^4 + (54*y^-6 - 351*y^-5 - 1561*y^-4 + 7398*y^-3 - 25007*y^-1 + 19467 + 19467*y - 25007*y^2 + 7398*y^4 - 1561*y^5 - 351*y^6 + 54*y^7)*q^5 + (394*y^-6 - 1561*y^-5 - 5916*y^-4 + 25007*y^-3 - 76461*y^-1 + 58537 + 58537*y - 76461*y^2 + 25007*y^4 - 5916*y^5 - 1561*y^6 + 394*y^7)*q^6 + O(q^7)]

AUTHORS:

- KENTA KOBAYASHI (2023-03-30): initial version

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

import math
from sage.all import LaurentPolynomialRing, LazyLaurentSeriesRing, QQ
from sage.WeakJacobiForm.weak_Jacobi_form.eisenstein import eisenstein

_R0 = LaurentPolynomialRing(QQ, "y")
_y = _R0.gen()

_R = LazyLaurentSeriesRing(_R0, "q")
_q = _R.gen()


def _function_t1(y):
    r"""

    Return a function which is proportional to the Jacobi theta function `\theta_1`.

    .. MATH::

    \begin{align}
    \theta_1(\tau, z)
    & = -\theta_{11}(\tau, z)
    = -i y^{1/2}q^{1/8}\sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2(n^2 + n)}\\
    t_1(y, q) & = \sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2(n^2 + n)}
    \end{align}
    where
    `q = e^{2 \pi i \tau}, y = e^{2 \pi i z}`

    INPUT:

    - ``y`` -- variable of the function which should be in ``LaurentPolynomialRing(QQ, "y")``.

    OUTPUT:

    the function `t_1` defined above.

    """

    # Return the pair of integral solutions of the equation "x^2 + x == n" or None.
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

    Return a function which is proportional to the Jacobi theta function `\theta_2`.

    .. MATH::

    \begin{align}
    \theta_2(\tau, z)
    & = \theta_{10}(\tau, z)
    = y^{1/2}q^{1/8}\sum_{n = -\infty}^\infty y^n q^{1/2(n^2 + n)} \\
    t_2(y, q) & = \sum_{n = -\infty}^\infty y^n q^{1/2(n^2 + n)}
    \end{align}

    INPUT:

    - ``y`` -- variable of the function which should be in ``LaurentPolynomialRing(QQ, "y")``.

    OUTPUT:

    the function `t_2` defined above.

    """

    # Return the pair of integral solutions of the equation "x^2 + x == n" or None.
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

    Return a function which is proportional to the Jacobi theta function `\theta_3`.

    .. MATH::

    \begin{align}
    \theta_3(\tau, z)
    & = \theta_{10}(\tau, z)
    = \sum_{n = -\infty}^\infty y^n q^{1/2n^2} \\
    t_3(y, q) & = \sum_{n = -\infty}^\infty y^n q^{n^2}
    \end{align}

    INPUT:

    - ``y`` -- variable of the function which should be in ``LaurentPolynomialRing(QQ, "y")``.

    OUTPUT:

    the function `t_3` defined above.

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

    Return a function which is proportional to the Jacobi theta function $`\theta_4`.

     .. MATH::

    \begin{align}
    \theta_4(\tau, z)
    & = \theta_{10}(\tau, z)
    = \sum_{n = -\infty}^\infty (-1)^n y^n q^{1/2n^2} \\
    t_4(y, q) & = \sum_{n = -\infty}^\infty (-1)^n y^n q^{n^2}
    \end{align}

    INPUT:

    - ``y`` -- variable of the function which should be in ``LaurentPolynomialRing(QQ, "y")``.

    OUTPUT:

    the function `t_4` defined above.

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

    Return a weak Jacobi form of weight `0` and index `1`.

    """
    f = 4 * (_function_t2(_y) / _function_t2(1)) ** 2 * _y

    g = (
        4 * (_function_t3(_y) / _function_t3(1)) ** 2
        + 4 * (_function_t4(_y) / _function_t4(1)) ** 2
    )
    return f + _R(lambda n: g.coefficient(2 * n), valuation=0)


def _phi_tilde_m2_1():
    r"""

    Return a weak Jacobi form of weight `-2` and index `1`.

    """
    h = _function_t3(1) * _function_t4(1)
    e3 = 1 / 2 * _function_t2(1) * _R(lambda n: h.coefficient(2 * n), valuation=0)
    return _function_t1(_y) ** 2 / e3**2 * _y


def _phi_tilde_0_3o2():
    r"""

    Return a weak Jacobi form of weight `0` and index `3/2` multiplied by `y^{1/2}`.

    """
    return _function_t1(_y**2) / _function_t1(_y) * _y


# nを2と3に分解するプログラム
def _decompose(n):
    r"""

    Return the list of all possible decompositions of ``n`` into sums of `2` and `3`.

    INPUT:

    - ``n`` -- an integer expected to be positive for meaningful calculations.

    OUTPUT:

    The set of pairs of coefficient of possible decompositions of ``n`` into sums of `2` and `3`,
    which is the set `\{(a, b) | 2a + 3b = n\}`.

    EXAMPLES::

        sage: from sage.WeakJacobiForm.weak_Jacobi_form.basis import _decompose
        sage: _decompose(1)
        set()
        sage: _decompose(4)
        {(2, 0)}
        sage: _decompose(8)
        {(1, 2), (4, 0)}

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

    Return a basis of the space of weak Jacobi forms of weight `0` and index ``index``.

    INPUT:

    - ``index`` -- integral index of weak Jacobi forms expected to be positive.

    OUTPUT:

    The list consisting of series in a basis of weight 0 and index ``index``.

    EXAMPLES::

        sage: from sage.WeakJacobiForm.weak_Jacobi_form.basis import basis_integral
        sage: basis_integral(1)
        [(y^-1 + 10 + y) + (10*y^-2 - 64*y^-1 + 108 - 64*y + 10*y^2)*q + (y^-3 + 108*y^-2 - 513*y^-1 + 808 - 513*y + 108*y^2 + y^3)*q^2 + (-64*y^-3 + 808*y^-2 - 2752*y^-1 + 4016 - 2752*y + 808*y^2 - 64*y^3)*q^3 + (10*y^-4 - 513*y^-3 + 4016*y^-2 - 11775*y^-1 + 16524 - 11775*y + 4016*y^2 - 513*y^3 + 10*y^4)*q^4 + (108*y^-4 - 2752*y^-3 + 16524*y^-2 - 43200*y^-1 + 58640 - 43200*y + 16524*y^2 - 2752*y^3 + 108*y^4)*q^5 + (y^-5 + 808*y^-4 - 11775*y^-3 + 58640*y^-2 - 141826*y^-1 + 188304 - 141826*y + 58640*y^2 - 11775*y^3 + 808*y^4 + y^5)*q^6 + O(q^7)]
        sage: basis_integral(3)
        [(y^-3 + 30*y^-2 + 303*y^-1 + 1060 + 303*y + 30*y^2 + y^3) + (30*y^-4 + 408*y^-3 - 456*y^-2 - 12696*y^-1 + 25428 - 12696*y - 456*y^2 + 408*y^3 + 30*y^4)*q + (303*y^-5 - 456*y^-4 - 14085*y^-3 + 143280*y^-2 - 477738*y^-1 + 697392 - 477738*y + 143280*y^2 - 14085*y^3 - 456*y^4 + 303*y^5)*q^2 + (1060*y^-6 - 12696*y^-5 + 143280*y^-4 - 1070296*y^-3 + 4380348*y^-2 - 9943440*y^-1 + 13003488 - 9943440*y + 4380348*y^2 - 1070296*y^3 + 143280*y^4 - 12696*y^5 + 1060*y^6)*q^3 + (303*y^-7 + 25428*y^-6 - 477738*y^-5 + 4380348*y^-4 - 22852404*y^-3 + 71497068*y^-2 - 139461585*y^-1 + 173777160 - 139461585*y + 71497068*y^2 - 22852404*y^3 + 4380348*y^4 - 477738*y^5 + 25428*y^6 + 303*y^7)*q^4 + (30*y^-8 - 12696*y^-7 + 697392*y^-6 - 9943440*y^-5 + 71497068*y^-4 - 304558128*y^-3 + 825547488*y^-2 - 1481082024*y^-1 + 1795708620 - 1481082024*y + 825547488*y^2 - 304558128*y^3 + 71497068*y^4 - 9943440*y^5 + 697392*y^6 - 12696*y^7 + 30*y^8)*q^5 + (y^-9 - 456*y^-8 - 477738*y^-7 + 13003488*y^-6 - 139461585*y^-5 + 825547488*y^-4 - 3056436603*y^-3 + 7532417376*y^-2 - 12788771979*y^-1 + 15228360016 - 12788771979*y + 7532417376*y^2 - 3056436603*y^3 + 825547488*y^4 - 139461585*y^5 + 13003488*y^6 - 477738*y^7 - 456*y^8 + y^9)*q^6 + O(q^7),
        (y^-3 + 6*y^-2 - 33*y^-1 + 52 - 33*y + 6*y^2 + y^3) + (6*y^-4 + 120*y^-3 + 2040*y^-2 - 9336*y^-1 + 14340 - 9336*y + 2040*y^2 + 120*y^3 + 6*y^4)*q + (-33*y^-5 + 2040*y^-4 - 30501*y^-3 + 169776*y^-2 - 436410*y^-1 + 590256 - 436410*y + 169776*y^2 - 30501*y^3 + 2040*y^4 - 33*y^5)*q^2 + (52*y^-6 - 9336*y^-5 + 169776*y^-4 - 1238200*y^-3 + 4575660*y^-2 - 9643728*y^-1 + 12291552 - 9643728*y + 4575660*y^2 - 1238200*y^3 + 169776*y^4 - 9336*y^5 + 52*y^6)*q^3 + (-33*y^-7 + 14340*y^-6 - 436410*y^-5 + 4575660*y^-4 - 23929812*y^-3 + 72597180*y^-2 - 137835345*y^-1 + 170028840 - 137835345*y + 72597180*y^2 - 23929812*y^3 + 4575660*y^4 - 436410*y^5 + 14340*y^6 - 33*y^7)*q^4 + (6*y^-8 - 9336*y^-7 + 590256*y^-6 - 9643728*y^-5 + 72597180*y^-4 - 310041072*y^-3 + 830703072*y^-2 - 1473696456*y^-1 + 1779000156 - 1473696456*y + 830703072*y^2 - 310041072*y^3 + 72597180*y^4 - 9643728*y^5 + 590256*y^6 - 9336*y^7 + 6*y^8)*q^5 + (y^-9 + 2040*y^-8 - 436410*y^-7 + 12291552*y^-6 - 137835345*y^-5 + 830703072*y^-4 - 3080138715*y^-3 + 7553501280*y^-2 - 12759364635*y^-1 + 15162554320 - 12759364635*y + 7553501280*y^2 - 3080138715*y^3 + 830703072*y^4 - 137835345*y^5 + 12291552*y^6 - 436410*y^7 + 2040*y^8 + y^9)*q^6 + O(q^7),
        (y^-3 - 6*y^-2 + 15*y^-1 - 20 + 15*y - 6*y^2 + y^3) + (-6*y^-4 - 456*y^-3 + 2856*y^-2 - 7224*y^-1 + 9660 - 7224*y + 2856*y^2 - 456*y^3 - 6*y^4)*q + (15*y^-5 + 2856*y^-4 - 40005*y^-3 + 182160*y^-2 - 414666*y^-1 + 539280 - 414666*y + 182160*y^2 - 40005*y^3 + 2856*y^4 + 15*y^5)*q^2 + (-20*y^-6 - 7224*y^-5 + 182160*y^-4 - 1325176*y^-3 + 4671156*y^-2 - 9491280*y^-1 + 11940768 - 9491280*y + 4671156*y^2 - 1325176*y^3 + 182160*y^4 - 7224*y^5 - 20*y^6)*q^3 + (15*y^-7 + 9660*y^-6 - 414666*y^-5 + 4671156*y^-4 - 24474996*y^-3 + 73142916*y^-2 - 137017041*y^-1 + 168165912 - 137017041*y + 73142916*y^2 - 24474996*y^3 + 4671156*y^4 - 414666*y^5 + 9660*y^6 + 15*y^7)*q^4 + (-6*y^-8 - 7224*y^-7 + 539280*y^-6 - 9491280*y^-5 + 73142916*y^-4 - 312795504*y^-3 + 833272224*y^-2 - 1469993736*y^-1 + 1770666660 - 1469993736*y + 833272224*y^2 - 312795504*y^3 + 73142916*y^4 - 9491280*y^5 + 539280*y^6 - 7224*y^7 - 6*y^8)*q^5 + (y^-9 + 2856*y^-8 - 414666*y^-7 + 11940768*y^-6 - 137017041*y^-5 + 833272224*y^-4 - 3092014395*y^-3 + 7564027680*y^-2 - 12744642603*y^-1 + 15129690352 - 12744642603*y + 7564027680*y^2 - 3092014395*y^3 + 833272224*y^4 - 137017041*y^5 + 11940768*y^6 - 414666*y^7 + 2856*y^8 + y^9)*q^6 + O(q^7)]

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

    Return a basis of the space of weak Jacobi forms of weight `0` and index ``double_index / 2`` for even ``double_index``,
    otherwise return a list of a basis multiplied by by `y^{1/2}`.

    INPUT:

    - `double_index` -- double of index of weak Jacobi forms expected to be greater than 1.

    OUTPUT:

    If `double_index` is even, the list consisting of series in a basis of weight `0` and index `double_index`.
    Otherwise, the list consisting of series in a basis of weight `0` and index `double_index` multiplied by $y^{1/2}$.

    EXAMPLE:

    For the case ``double_index`` is even, the output is equal to ``basis_integral(int(double_index / 2))``

        sage: from sage.WeakJacobiForm.weak_Jacobi_form.basis import basis_half_integral
        sage: basis_half_integral(4)
        [(y^-2 + 20*y^-1 + 102 + 20*y + y^2) + (20*y^-3 + 72*y^-2 - 1044*y^-1 + 1904 - 1044*y + 72*y^2 + 20*y^3)*q + (102*y^-4 - 1044*y^-3 + 7392*y^-2 - 23532*y^-1 + 34164 - 23532*y + 7392*y^2 - 1044*y^3 + 102*y^4)*q^2 + (20*y^-5 + 1904*y^-4 - 23532*y^-3 + 115552*y^-2 - 283688*y^-1 + 379488 - 283688*y + 115552*y^2 - 23532*y^3 + 1904*y^4 + 20*y^5)*q^3 + (y^-6 - 1044*y^-5 + 34164*y^-4 - 283688*y^-3 + 1107711*y^-2 - 2402244*y^-1 + 3090200 - 2402244*y + 1107711*y^2 - 283688*y^3 + 34164*y^4 - 1044*y^5 + y^6)*q^4 + (72*y^-6 - 23532*y^-5 + 379488*y^-4 - 2402244*y^-3 + 8066472*y^-2 - 16135248*y^-1 + 20229984 - 16135248*y + 8066472*y^2 - 2402244*y^3 + 379488*y^4 - 23532*y^5 + 72*y^6)*q^5 + (20*y^-7 + 7392*y^-6 - 283688*y^-5 + 3090200*y^-4 - 16135248*y^-3 + 48552352*y^-2 - 91617180*y^-1 + 112772304 - 91617180*y + 48552352*y^2 - 16135248*y^3 + 3090200*y^4 - 283688*y^5 + 7392*y^6 + 20*y^7)*q^6 + O(q^7),
        (y^-2 - 4*y^-1 + 6 - 4*y + y^2) + (-4*y^-3 + 264*y^-2 - 1020*y^-1 + 1520 - 1020*y + 264*y^2 - 4*y^3)*q + (6*y^-4 - 1020*y^-3 + 8160*y^-2 - 23556*y^-1 + 32820 - 23556*y + 8160*y^2 - 1020*y^3 + 6*y^4)*q^2 + (-4*y^-5 + 1520*y^-4 - 23556*y^-3 + 117856*y^-2 - 283640*y^-1 + 375648 - 283640*y + 117856*y^2 - 23556*y^3 + 1520*y^4 - 4*y^5)*q^3 + (y^-6 - 1020*y^-5 + 32820*y^-4 - 283640*y^-3 + 1113855*y^-2 - 2402316*y^-1 + 3080600 - 2402316*y + 1113855*y^2 - 283640*y^3 + 32820*y^4 - 1020*y^5 + y^6)*q^4 + (264*y^-6 - 23556*y^-5 + 375648*y^-4 - 2402316*y^-3 + 8081256*y^-2 - 16135152*y^-1 + 20207712 - 16135152*y + 8081256*y^2 - 2402316*y^3 + 375648*y^4 - 23556*y^5 + 264*y^6)*q^5 + (-4*y^-7 + 8160*y^-6 - 283640*y^-5 + 3080600*y^-4 - 16135152*y^-3 + 48585376*y^-2 - 91617300*y^-1 + 112723920 - 91617300*y + 48585376*y^2 - 16135152*y^3 + 3080600*y^4 - 283640*y^5 + 8160*y^6 - 4*y^7)*q^6 + O(q^7)]

    For the case ``double_index`` is odd, the output multiplied by `y^{-1/2}` is the correct basis.

        sage: basis_half_integral(3)
        [(1 + y) + (-y^-2 + 1 + y - y^3)*q + (-y^-3 - y^-2 + 2 + 2*y - y^3 - y^4)*q^2 + (-y^-3 - 2*y^-2 + 3 + 3*y - 2*y^3 - y^4)*q^3 + (-2*y^-3 - 3*y^-2 + 5 + 5*y - 3*y^3 - 2*y^4)*q^4 + (y^-5 - 3*y^-3 - 5*y^-2 + 7 + 7*y - 5*y^3 - 3*y^4 + y^6)*q^5 + (y^-5 - 5*y^-3 - 7*y^-2 + 11 + 11*y - 7*y^3 - 5*y^4 + y^6)*q^6 + O(q^7)]
        sage: basis_half_integral(5)
        [(y^-1 + 11 + 11*y + y^2) + (-y^-3 - 54*y^-1 + 55 + 55*y - 54*y^2 - y^4)*q + (-11*y^-4 + 54*y^-3 - 394*y^-1 + 351 + 351*y - 394*y^2 + 54*y^4 - 11*y^5)*q^2 + (-11*y^-5 - 55*y^-4 + 394*y^-3 - 1889*y^-1 + 1561 + 1561*y - 1889*y^2 + 394*y^4 - 55*y^5 - 11*y^6)*q^3 + (-y^-6 - 55*y^-5 - 351*y^-4 + 1889*y^-3 - 7398*y^-1 + 5916 + 5916*y - 7398*y^2 + 1889*y^4 - 351*y^5 - 55*y^6 - y^7)*q^4 + (54*y^-6 - 351*y^-5 - 1561*y^-4 + 7398*y^-3 - 25007*y^-1 + 19467 + 19467*y - 25007*y^2 + 7398*y^4 - 1561*y^5 - 351*y^6 + 54*y^7)*q^5 + (394*y^-6 - 1561*y^-5 - 5916*y^-4 + 25007*y^-3 - 76461*y^-1 + 58537 + 58537*y - 76461*y^2 + 25007*y^4 - 5916*y^5 - 1561*y^6 + 394*y^7)*q^6 + O(q^7)]

        For the case ``double_index`` is `1`, the output is empty.

        sage: basis_half_integral(1)
        []

    """
    if double_index == 1:
        return []
    elif double_index % 2 == 0:
        return basis_integral(int(double_index / 2))
    else:
        index = int((double_index - 3) / 2)
        return list(map(lambda f: f * _phi_tilde_0_3o2(), basis_integral(index)))
