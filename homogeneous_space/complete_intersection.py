r"""
Class of varieties which is the zeros of a general section of an equivariant vector bundle of a homogeneous space
================================================

This module contains a class:
    - ``CompleteIntersection`` -- an implementation of ``IVariety``, representing the zeros of a general section of an equivariant vector bundle of a homogeneous space
This is specialized in computing Chern characters and Todd classes from Chern classes.

EXAMPLE:
    sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
    sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
    sage: from sage.EllipticGenus.homogeneous_space.complete_intersection import CompleteIntersection
    sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
    sage: X = HomogeneousSpace(P)
    sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
    sage: quintic = CompleteIntersection(X, E)
    sage: quintic.numerical_integration_by_localization(quintic.chern_classes()[3])
    -200


AUTHORS:

- KENTA KOBAYASHI (2023-04-04): initial version

REFERENCES:

.. [IMOU2022] Atsushi Ito, Makoto Miura, Shinnosuke Okawa, and Kazushi Ueda. Calabi--yau complete intersec-
tions in exceptional grassmannians, 2022.


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

from sage.all import prod, vector
from sage.EllipticGenus.homogeneous_space.homogeneous_space import (
    EquivariantVectorBundle,
    HomogeneousSpace,
)
from sage.EllipticGenus.homogeneous_space.interfaces import IVariety, IVectorBundle

homogeneous_part = lambda F, degree: sum(
    c * m for c, m in F if m.total_degree() == degree
)


class CompleteIntersection(IVariety):
    r"""

    Class representing the zeros of a general section of an equivariant vector bundle of a homogeneous space

    """

    def __init__(
        self,
        homogeneous_space: HomogeneousSpace,
        vector_bundle: EquivariantVectorBundle,
    ) -> None:
        r"""

        Constructor of this class

        INPUT:

        - ``homogeneous_space`` -- ``HomogeneousSpace`` -- the variety contains this

        - ``vector_bundle`` -- ``EquivariantVectorBundle``


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
            sage: from sage.EllipticGenus.homogeneous_space.complete_intersection import CompleteIntersection
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
            sage: CompleteIntersection(X, E)
            a complete intersection of a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1] and an equivariant vector bundle on a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1] associated to (5, 0, 0, 0, 0)

        """
        self.homogeneous_space = homogeneous_space
        self.vector_bundle = vector_bundle
        self.dim = self.homogeneous_space.dimension() - self.vector_bundle.rank()

    def __repr__(self) -> str:
        return f"a complete intersection of {self.homogeneous_space} and {self.vector_bundle}"

    def dimension(self) -> int:
        r"""
        Return the dimension of this variety

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
            sage: from sage.EllipticGenus.homogeneous_space.complete_intersection import CompleteIntersection
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
            sage: quintic = CompleteIntersection(X, E)
            sage: quintic.dimension()
            3
        """
        return self.dim

    def tangent_bundle(self):
        r"""
        Return the tangent bundle of this variety

        """
        ci = self

        class VB(IVectorBundle):
            def base(self) -> IVariety:
                return ci

            def rank(self) -> int:
                return ci.dim

            def chern_classes(self) -> list:
                return ci.chern_classes()

        return VB()

    def chern_classes(self):
        r"""
        Return the list of homogeneous parts of Chern classes of the tangent bundle of this variety

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
            sage: from sage.EllipticGenus.homogeneous_space.complete_intersection import CompleteIntersection
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
            sage: quintic = CompleteIntersection(X, E)
            sage: quintic.chern_classes()
            [1,
             -x0 - x1 - x2 - x3 - x4,
             11*x0^2 + 2*x0*x1 + 2*x0*x2 + x1*x2 + 2*x0*x3 + x1*x3 + x2*x3 + 2*x0*x4 + x1*x4 + x2*x4 + x3*x4,
             -51*x0^3 - 13*x0^2*x1 - 13*x0^2*x2 - 3*x0*x1*x2 - 13*x0^2*x3 - 3*x0*x1*x3 - 3*x0*x2*x3 - x1*x2*x3 - 13*x0^2*x4 - 3*x0*x1*x4 - 3*x0*x2*x4 - x1*x2*x4 - 3*x0*x3*x4 - x1*x3*x4 - x2*x3*x4]
        """

        def class_from_weight(weight):
            return sum(
                weight[i] * self.homogeneous_space.x[i]
                for i in range(
                    self.homogeneous_space.parabolic_subgroup.ambient_space().dimension()
                )
            )

        def geometric_sequence(n, x):
            return sum(x**i for i in range(0, n + 1))

        cc = prod(1 + x for x in self.homogeneous_space.tangent_weights) * prod(
            geometric_sequence(self.dim, -class_from_weight(vector(w))) ** i
            for w, i in self.vector_bundle.weight_multiplicities.items()
        )

        return [homogeneous_part(cc, i) for i in range(0, self.dim + 1)]

    def numerical_integration_by_localization(self, f):
        r"""

        Return the numerical computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the numerical computation of the integration of the  equivariant cohomology class ``f``.


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
            sage: from sage.EllipticGenus.homogeneous_space.complete_intersection import CompleteIntersection
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(5, 0, 0, 0, 0))
            sage: quintic = CompleteIntersection(X, E)
            sage: quintic.numerical_integration_by_localization(quintic.chern_classes()[3])
            -200

        """
        top_of_f = homogeneous_part(f, self.dim)
        c_top = (
            self.vector_bundle.chern_classes()[self.vector_bundle.rank()]
            if self.dim >= 0
            else 0
        )
        return self.homogeneous_space.numerical_integration_by_localization(
            top_of_f * c_top
        )

    def integration(self, f) -> int:
        r"""

        Implementation of the abstract method.

        """
        return self.numerical_integration_by_localization(f)
