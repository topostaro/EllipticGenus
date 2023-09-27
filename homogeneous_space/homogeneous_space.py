r"""
Classes of homogeneous spaces and equivariant vector bundles
============================================================

This module contains classes:

- ``HomogeneousSpace`` -- a class of homogeneous spaces inherit ``AlmostComplexManifold``,
- ``EquivariantVectorBundle`` -- a class of equivariant vector bundles on homogeneous spaces, which inherit ``VectorBundle``.
- ``IrreducibleEquivariantVectorBundle`` -- a class of irreducible equivariant vector bundles on homogeneous spaces, which inherit ``EquivariantVectorBundle``.

These are specialized in computing Chern characters and Todd classes from Chern classes.

EXAMPLES::

    sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, EquivariantVectorBundle
    sage: from homogeneous_space.parabolic import ParabolicSubgroup
    sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
    sage: X = HomogeneousSpace(P)
    sage: X
    a homogeneous_space associated to the parabolic subgroup of ['A', 4] with crossed-out nodes [1]
    sage: E = EquivariantVectorBundle(X,{(2, 0, 0, 0, 0): 1, (3, 0, 0, 0, 0): 1})
    sage: E.chern_classes()
    [1, 5*x0, 6*x0^2, 0, 0]


AUTHORS:

- KENTA KOBAYASHI, AKIHITO NAKAMURA and KAZUSHI UEDA (2023-04-04): initial version


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
from sage.all import PolynomialRing, QQ, prod, RealField, vector, WeylGroup, Matrix
from homogeneous_space.parabolic import ParabolicSubgroup
from homogeneous_space.interfaces import AlmostComplexManifold, VectorBundle
from functools import cache


# `degree`次部分を取り出す関数
homogeneous_part = lambda F, degree: sum(
    c * m for c, m in F if m.total_degree() == degree
)


class HomogeneousSpace(AlmostComplexManifold):
    r"""

    Class representing a homogeneous space, which is the quotient by a parabolic subgroup.

    """

    def __init__(self, parabolic_subgroup: ParabolicSubgroup) -> None:
        r"""

        Constructor of this class

        INPUT:

        - ``parabolic_subgroup`` -- ``ParabolicSubgroup``

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X.dimension()
            4

        """
        return self.dim

    def tangent_bundle(self) -> VectorBundle:
        r"""
        Return the tangent bundle of this variety
        """

        tangent_weights = {
            r: 1
            for r in set(self.parabolic_subgroup.R_G.positive_roots())
            - set(self.parabolic_subgroup.positive_roots())
        }

        return EquivariantVectorBundle(self, tangent_weights)

    # 次数ごとのchern類
    @cache
    def chern_classes(self):
        r"""
        Return the list of homogeneous parts of Chern classes of the tangent bundle of this variety

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

    @cache
    def numerical_integration_by_localization(self, f):
        r"""

        Return the numerical computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the numerical computation of the integration of the  equivariant cohomology class ``f``.

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

        len_of_wg_of_L = (
            1
            if self.parabolic_subgroup.L == None
            else len(WeylGroup(self.parabolic_subgroup.L))
        )

        if top_of_f == 0:
            return 0
        else:
            return (
                sum(
                    [
                        top_of_f(x) / denominator_in_localization(x)
                        for x in orbit_of_random_x
                    ]
                ).round()
                / len_of_wg_of_L
            )

    @cache
    def symbolic_integration_by_localization(self, f):
        r"""

        Return the symbolic computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the symbolic computation of the integration of the  equivariant cohomology class ``f``.

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: X.symbolic_integration_by_localization(X.chern_classes()[X.dimension()])
            5

        """
        S_G = WeylGroup(
            self.parabolic_subgroup.G.root_system(), implementation="permutation"
        )

        s_G = S_G.simple_reflections()
        s_P = [
            s_G[i]
            for i in set(range(1, len(s_G) + 1))
            - set(self.parabolic_subgroup.crossed_out_nodes)
        ]
        S_P = S_G.subgroup(s_P)

        def cosetRep(G, H):
            rep = []
            g = set(G.list())
            h = set(G(e) for e in H.list())
            coset = set()
            while g:
                p = g.pop()
                rep.append(p)
                coset = set(p * e for e in h)
                g = g - coset
                coset.clear()
                if len(rep) * H.order() == G.order():
                    break

            return rep

        W_G = WeylGroup(self.parabolic_subgroup.G.root_system())
        W_G_mod_W_P = [W_G(w) for w in cosetRep(S_G, S_P)]

        def weyl_group_action(w, f):
            return f(*(w.inverse() * vector(self.ring, self.ring.gens())))

        return sum(
            [
                weyl_group_action(w, f)
                / weyl_group_action(w, prod(self.tangent_weights))
                for w in W_G_mod_W_P
            ]
        )

    def integration(self, f, option="symbolic") -> int:
        r"""

        Implementation of the abstract method.

        """

        if option == "symbolic":
            return self.symbolic_integration_by_localization(f)
        if option == "numerical":
            return self.numerical_integration_by_localization(f)
        else:
            raise TypeError(
                f"Invalid option in integration on HomogeneousSpace: {option}"
            )


class EquivariantVectorBundle(VectorBundle):
    r"""

    Class representing an equivariant vector bundle on a homogeneous space

    This class contains the base homogeneous space `G/P` and a representation of `G`.

    """

    def __init__(
        self, homogeneous_space: HomogeneousSpace, weight_multiplicities: dict
    ) -> None:
        r"""

        Constructor of this class

        This constructor takes the base homogeneous space `G/P` and weights of `G` with their multiplicities. This construct an equivariant vector bundle associated to the representation of `G` associated to the weights with multiplicities.

        INPUT:

        - ``homogeneous_space`` -- ``HomogeneousSpace`` -- the base space of this vector bundle

        - ``weight_multiplicities`` -- dictionary from weights to their multiplicities

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, EquivariantVectorBundle
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, EquivariantVectorBundle
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, EquivariantVectorBundle
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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

        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, EquivariantVectorBundle
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
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


class CompletelyReducibleEquivariantVectorBundle(EquivariantVectorBundle):
    def __init__(
        self, homogeneous_space: HomogeneousSpace, highest_weights: list
    ) -> None:
        if len(highest_weights) == 1:
            self.is_irr = True
        else:
            self.is_irr = False

        def union_multisets(mset1, mset2):
            result = dict(mset1)

            for key, value in mset2.items():
                if key in result:
                    result[key] += value
                else:
                    result[key] = value

            return result

        def flatten(list_of_mset):
            result = dict(list_of_mset[0])

            for mset in list_of_mset[1:]:
                result = union_multisets(result, mset)

            return result

        weight_multiplicities = flatten(
            list(
                map(
                    lambda w: homogeneous_space.parabolic_subgroup.weight_multiplicities(
                        w
                    ),
                    highest_weights,
                )
            )
        )

        self.highest_weights = highest_weights
        super().__init__(homogeneous_space, weight_multiplicities)

    def is_irreducible(self) -> bool:
        return self.is_irr


class IrreducibleEquivariantVectorBundle(CompletelyReducibleEquivariantVectorBundle):
    def __init__(self, homogeneous_space, weight) -> None:
        r"""

        Constructor of this class

        This constructor takes the base homogeneous space `G/P` and weights of `G` with their multiplicities. This construct an equivariant vector bundle associated to the highest weight representation of `G` associated to the weight ``weight``.

        INPUT:

        - ``weight`` -- list of integers -- this represents the coefficients of the fundamental weights


        EXAMPLES::

            sage: from homogeneous_space.homogeneous_space import HomogeneousSpace, IrreducibleEquivariantVectorBundle
            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
            sage: X = HomogeneousSpace(P)
            sage: E = IrreducibleEquivariantVectorBundle(X,(0, 1, 0, 0, 0))
            sage: E.rank()
            4

        """
        self.weight = weight
        super().__init__(homogeneous_space, [weight])

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.homogeneous_space} associated to {self.weight}"

    def is_irreducible(self) -> bool:
        return True
