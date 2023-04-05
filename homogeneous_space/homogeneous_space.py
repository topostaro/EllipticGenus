r"""
Classes of homogeneous spaces and equivariant vector bundles
================================================

This module contains classes:
    - ``HomogeneousSpace`` -- a class of homogeneous spaces inherit ``IVariety``,
    - ``EquivariantVectorBundle`` -- a class of equivariant vector bundles on homogeneous spaces, which inherit ``IVectorBundle``.
    - ``IrreducibleEquivariantVectorBundle`` -- a class of irreducible equivariant vector bundles on homogeneous spaces, which inherit ``EquivariantVectorBundle``.
These are specialized in computing Chern characters and Todd classes from Chern classes.

EXAMPLE:
    sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
    sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
    sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
    sage: X = HomogeneousSpace(P)
    sage: X
    a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1]
    sage: E = EquivariantVectorBundle(X,{(2, 0, 0, 0, 0): 1, (3, 0, 0, 0, 0): 1})
    sage: E.chern_classes()
    [1, 5*x0, 6*x0^2, 0, 0]


AUTHORS:

- KENTA KOBAYASHI (2023-04-04): initial version

.. [BE1989] Robert J. Baston and Michael G. Eastwood, The Penrose transform, Oxford Mathematical Monographs,
The Clarendon Press, Oxford University Press, New York, 1989, Its interaction with representation
theory, Oxford Science Publications.

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

from random import random
from sage.all import PolynomialRing, QQ, prod, RealField, vector, WeylGroup
from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
from sage.EllipticGenus.homogeneous_space.interfaces import IVariety, IVectorBundle

# `degree`次部分を取り出す関数
homogeneous_part = lambda F, degree: sum(
    c * m for c, m in F if m.total_degree() == degree
)


class HomogeneousSpace(IVariety):
    r"""

    Class representing a homogeneous space, which is the quotient by a parabolic subgroup.

    """

    def __init__(self, parabolic_subgroup: ParabolicSubgroup) -> None:
        r"""

        Constructor of this class

        INPUT:

        - ``parabolic_subgroup`` -- ``ParabolicSubgroup``


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X
            a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1]

        """
        self.parabolic_subgroup = parabolic_subgroup

        # コホモロジー環を含む環
        self.ring = PolynomialRing(
            QQ, "x", parabolic_subgroup.ambient_space().dimension()
        )
        self.x = self.ring.gens()

        self.tangent_weights = [
            sum(
                r[l] * self.x[l]
                for l in range(parabolic_subgroup.ambient_space().dimension())
            )
            for r in set(parabolic_subgroup.R_G.positive_roots())
            - set(parabolic_subgroup.positive_roots())
        ]
        self.dim = len(self.tangent_weights)

    def __repr__(self) -> str:
        return f"a homogeneous_space associated to {self.parabolic_subgroup}"

    def dimension(self) -> int:
        r"""
        Return the dimension of this variety

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X.dimension()
            4

        """
        return self.dim

    def tangent_bundle(self) -> IVectorBundle:
        r"""
        Return the tangent bundle of this variety
        """
        tangent_weights = {
            w: 1
            for w in set(self.parabolic_subgroup.R_G.positive_roots())
            - set(self.parabolic_subgroup.positive_roots())
        }
        return EquivariantVectorBundle(self, tangent_weights)

    # 次数ごとのchern類
    def chern_classes(self):
        r"""
        Return the list of homogeneous parts of Chern classes of the tangent bundle of this variety


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X.chern_classes()
            [1,
             4*x0 - x1 - x2 - x3 - x4,
             6*x0^2 - 3*x0*x1 - 3*x0*x2 + x1*x2 - 3*x0*x3 + x1*x3 + x2*x3 - 3*x0*x4 + x1*x4 + x2*x4 + x3*x4,
             4*x0^3 - 3*x0^2*x1 - 3*x0^2*x2 + 2*x0*x1*x2 - 3*x0^2*x3 + 2*x0*x1*x3 + 2*x0*x2*x3 - x1*x2*x3 - 3*x0^2*x4 + 2*x0*x1*x4 + 2*x0*x2*x4 - x1*x2*x4 + 2*x0*x3*x4 - x1*x3*x4 - x2*x3*x4,
             x0^4 - x0^3*x1 - x0^3*x2 + x0^2*x1*x2 - x0^3*x3 + x0^2*x1*x3 + x0^2*x2*x3 - x0*x1*x2*x3 - x0^3*x4 + x0^2*x1*x4 + x0^2*x2*x4 - x0*x1*x2*x4 + x0^2*x3*x4 - x0*x1*x3*x4 - x0*x2*x3*x4 + x1*x2*x3*x4]

        """
        return [
            homogeneous_part(prod(1 + x for x in self.tangent_weights), i)
            for i in (range(0, self.dim + 1))
        ]

    def numerical_integration_by_localization(self, f):
        r"""

        Return the numerical computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the numerical computation of the integration of the  equivariant cohomology class ``f``.

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X.numerical_integration_by_localization(X.chern_classes()[X.dimension()])
            5

        """

        random_x = [
            RealField(1000)(random())
            for i in range(self.parabolic_subgroup.ambient_space().dimension())
        ]
        orbit_of_random_x = [
            (w.inverse() * vector(RealField(1000), random_x)).list()
            for w in WeylGroup(self.parabolic_subgroup.G)
        ]
        top_of_f = homogeneous_part(f, self.dim)
        denominator_in_localization = prod(self.tangent_weights)
        if top_of_f == 0:
            return 0
        else:
            return sum(
                [
                    top_of_f(x) / denominator_in_localization(x)
                    for x in orbit_of_random_x
                ]
            ).round() / len(WeylGroup(self.parabolic_subgroup.L))

    def integration(self, f) -> int:
        r"""

        Implementation of the abstract method.

        """

        return self.numerical_integration_by_localization(f)


class EquivariantVectorBundle(IVectorBundle):
    r"""

    Class representing a equivariant vector bundle on a homogeneous space

    This class contains the base homogeneous space `G/P` and a representation of `G`.

    """

    def __init__(self, homogeneous_space, weight_multiplicities) -> None:
        r"""

        Constructor of this class

        This constructor takes the base homogeneous space `G/P` and weights of `G` with their multiplicities. This construct an equivariant vector bundle associated to the representation of `G` associated to the weights with multiplicities.

        INPUT:

        - ``homogeneous_space`` -- ``HomogeneousSpace`` -- the base space of this vector bundle

        - ``weight_multiplicities`` -- dictionary from weights to their multiplicities

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: EquivariantVectorBundle(X,{(3, 2, 0, 0, 0): 1, (3, 1, 1, 0, 0): 1, (3, 0, 2, 0, 0): 1, (3, 1, 0, 1, 0): 1, (3, 0, 1, 1, 0): 1, (3, 1, 0, 0, 1): 1, (3, 0, 0, 2, 0): 1, (3, 0, 1, 0, 1): 1, (3, 0, 0, 1, 1): 1, (3, 0, 0, 0, 2): 1, (4, 0, 0, 0, 0): 1,})
            an equivariant vector bundle on a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1] associated to {(3, 2, 0, 0, 0): 1, (3, 1, 1, 0, 0): 1, (3, 0, 2, 0, 0): 1, (3, 1, 0, 1, 0): 1, (3, 0, 1, 1, 0): 1, (3, 1, 0, 0, 1): 1, (3, 0, 0, 2, 0): 1, (3, 0, 1, 0, 1): 1, (3, 0, 0, 1, 1): 1, (3, 0, 0, 0, 2): 1, (4, 0, 0, 0, 0): 1}

        """
        self.homogeneous_space = homogeneous_space
        self.weight_multiplicities = weight_multiplicities
        self.rk = sum(v for v in self.weight_multiplicities.values())

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.homogeneous_space} associated to {self.weight_multiplicities}"

    def rank(self) -> int:
        r"""
        Return the rank of this vector bundle


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = EquivariantVectorBundle(X,{(3, 2, 0, 0, 0): 1, (3, 1, 1, 0, 0): 1, (3, 0, 2, 0, 0): 1, (3, 1, 0, 1, 0): 1, (3, 0, 1, 1, 0): 1, (3, 1, 0, 0, 1): 1, (3, 0, 0, 2, 0): 1, (3, 0, 1, 0, 1): 1, (3, 0, 0, 1, 1): 1, (3, 0, 0, 0, 2): 1, (4, 0, 0, 0, 0): 1,})
            sage: E.rank()
            11

        """
        return self.rk

    def base(self) -> HomogeneousSpace:
        r"""
        Return the base space of this vector bundle

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = EquivariantVectorBundle(X,{(3, 2, 0, 0, 0): 1, (3, 1, 1, 0, 0): 1, (3, 0, 2, 0, 0): 1, (3, 1, 0, 1, 0): 1, (3, 0, 1, 1, 0): 1, (3, 1, 0, 0, 1): 1, (3, 0, 0, 2, 0): 1, (3, 0, 1, 0, 1): 1, (3, 0, 0, 1, 1): 1, (3, 0, 0, 0, 2): 1, (4, 0, 0, 0, 0): 1,})
            sage: E.base()
            a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1]
        """
        return self.homogeneous_space

    def chern_classes(self):
        r"""
        Return the list of homogeneous parts of Chern classes of this vector bundle

        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = EquivariantVectorBundle(X,{(2, 0, 0, 0, 0): 1, (3, 0, 0, 0, 0): 1})
            sage: E.chern_classes()
            [1, 5*x0, 6*x0^2, 0, 0]

        """

        def class_from_weight(weight):
            return sum(
                weight[i] * self.homogeneous_space.x[i]
                for i in range(
                    self.homogeneous_space.parabolic_subgroup.ambient_space().dimension()
                )
            )

        cc = prod(
            (1 + class_from_weight(vector(w))) ** i
            for w, i in self.weight_multiplicities.items()
        )

        return [
            homogeneous_part(cc, i) for i in (range(0, self.homogeneous_space.dim + 1))
        ]


class IrreducibleEquivariantVectorBundle(EquivariantVectorBundle):
    def __init__(self, homogeneous_space, weight) -> None:
        r"""

        Constructor of this class

        This constructor takes the base homogeneous space `G/P` and weights of `G` with their multiplicities. This construct an equivariant vector bundle associated to the highest weight representation of `G` associated to the weight ``weight``.

        INPUT:

        - ``weight`` -- list of integers -- this represents the coefficients of the fundamental weights


        EXAMPLE:
            sage: from sage.EllipticGenus.homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from sage.EllipticGenus.homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(0, 1, 0, 0, 0))
            sage: E.rank()
            4

        """
        self.weight = weight
        super().__init__(
            homogeneous_space,
            homogeneous_space.parabolic_subgroup.weight_multiplicities(weight),
        )

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.homogeneous_space} associated to {self.weight}"
