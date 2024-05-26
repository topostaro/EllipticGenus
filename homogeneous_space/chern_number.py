r"""
Computation of Chern numbers of varieties
================================================

This module implements a computation of Chern numbers of varieties.

EXAMPLES::

    sage: from homogeneous_space.parabolic import ParabolicSubgroup
    sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
    sage: from homogeneous_space.complete_intersection import CompleteIntersection
    sage: from homogeneous_space.chern_number import chern_number
    sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
    sage: X = HomogeneousSpace(P)
    sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
    sage: quintic = CompleteIntersection(E)
    sage: chern_number(quintic, [1, 1, 1])
    0
    sage: chern_number(quintic, [3])
    -200


AUTHORS:

- KENTA KOBAYASHI (2023-04-04): initial version

REFERENCES:


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

from sage.all import prod
from homogeneous_space.interfaces import AlmostComplexManifold


# `degree`次部分を取り出す関数
def homogeneous_part(F, degree):
    return sum(
        c * m for c, m in F if m.total_degree() == degree
    )


def chern_number(
    manifold: AlmostComplexManifold, degrees: list, option="symbolic"
) -> int:
    r"""

    Return the Chern number of ``manifold`` with degree ``degree``
    by integration with an option ``option``

    INPUT:

    - ``manifold`` -- ``AlmostComplexManifold``

    - ``degree`` -- list of integers -- this list represents the integrant. For example, ``[1,1,3]`` represents `c_1 c_1 c_3`.

    - ``option`` -- string -- this specifying the integration option.

    OUTPUT:

    the Chern number of ``manifold`` with degree ``degree``

    EXAMPLES::

        sage: from homogeneous_space.parabolic import ParabolicSubgroup
        sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
        sage: from homogeneous_space.complete_intersection import CompleteIntersection
        sage: from homogeneous_space.chern_number import chern_number
        sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
        sage: X = HomogeneousSpace(P)
        sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
        sage: quintic = CompleteIntersection(E)
        sage: chern_number(quintic, [1, 1, 1])
        0
        sage: chern_number(quintic, [3])
        -200

    """
    if sum(d for d in degrees) != manifold.dimension():
        return 0
    else:
        chern_classes = manifold.chern_classes()
        return manifold.integration(prod([chern_classes[d] for d in degrees]), option)
