r"""
Computation of Euler characteristic of complex vector bundles
=============================================================

This module implements a computation of Euler characteristic of complex vector bundles.

EXAMPLES::

    sage: from homogeneous_space.parabolic import ParabolicSubgroup
    sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
    sage: from homogeneous_space.euler_characteristic import euler_characteristic
    sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
    sage: X = HomogeneousSpace(P)
    sage: E = IrreducibleEquivariantVectorBundle(X,(3, 0, 0, 0, 0))
    sage: euler_characteristic(X, E)
    35


AUTHORS:

- KENTA KOBAYASHI (2023-04-04): initial version

REFERENCES:

.. [Hir1978] Friedrich Hirzebruch, Topological methods in algebraic geometry, Springer (1978)


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

from homogeneous_space.interfaces import AlmostComplexManifold, VectorBundle


def euler_characteristic(
    manifold: AlmostComplexManifold, vector_bundle: VectorBundle
) -> int:
    r"""

    Return the Euler characteristic of ``vector_bundle`` on ``manifold``

    INPUT:

    - ``manifold`` -- object of ``AlmostComplexManifold`` -- the base space.

    - ``vector_bundle`` -- object of ``VectorBundle``

    OUTPUT:

    the Euler characteristic of ``vector_bundle`` on ``manifold``, that is computed as the integration of the multiplication of the Chern character of ``vector_bundle`` and the Todd genus of ``manifold``.

    ..MATH::

        \chi (X, E) = \int_X ch(E) td(X)

    EXAMPLES::

        sage: from homogeneous_space.parabolic import ParabolicSubgroup
        sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
        sage: from homogeneous_space.euler_characteristic import euler_characteristic
        sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
        sage: X = HomogeneousSpace(P)
        sage: E = IrreducibleEquivariantVectorBundle(X,(3, 0, 0, 0, 0))
        sage: euler_characteristic(X, E)
        35

    """
    chern_character = vector_bundle.chern_character()
    todd_classes = manifold.todd_classes()

    return manifold.integration(
        sum(
            chern_character[i] * todd_classes[manifold.dimension() - i]
            for i in range(manifold.dimension() + 1)
        )
    )
