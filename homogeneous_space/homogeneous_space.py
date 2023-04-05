r"""
Classes of homogeneous spaces and equivariant vector bundles
================================================

This module contains classes:
    - ``HomogeneousSpace`` -- a class of homogeneous spaces inherit ``IVariety``,
    - ``EquivariantVectorBundle`` -- a class of equivariant vector bundles on homogeneous spaces, which inherit ``IVectorBundle``.
    - ``IrreducibleEquivariantVectorBundle`` -- a class of irreducible equivariant vector bundles on homogeneous spaces, which inherit ``EquivariantVectorBundle``.
These are specialized in computing Chern characters and Todd classes from Chern classes.

EXAMPLE:


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

from sage.all import PolynomialRing
from parabolic import ParabolicSubgroup
from interfaces import *


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

        """
        self.parabolic_subgroup = parabolic_subgroup

        # コホモロジー環を含む環
        self.ring = PolynomialRing(
            QQ, "x", parabolic_subgroup.ambient_space_dimension()
        )
        self.x = self.ring.gens()

        self.tangent_weights = [
            sum(
                r[l] * self.x[l]
                for l in range(parabolic_subgroup.ambient_space_dimension())
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
        """
        return [
            homogeneous_part(prod(1 + x for x in self.tangent_weights), i)
            for i in ((0.0).self.dim)
        ]

    def numerical_integration_by_localization(self, f):
        r"""

        Return the numerical computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the numerical computation of the integration of the  equivariant cohomology class ``f``.

        """
        random_x = [
            RealField(1000)(random())
            for i in range(self.parabolic_subgroup.ambient_space_dimension())
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

        """
        self.homogeneous_space = homogeneous_space
        self.weight_multiplicities = weight_multiplicities
        self.rk = sum(v for v in self.weight_multiplicities.values())

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.homogeneous_space} associated to {self.weight_multiplicities}"

    def rank(self) -> int:
        r"""
        Return the rank of this vector bundle
        """
        return self.rk

    def base(self) -> HomogeneousSpace:
        r"""
        Return the base space of this vector bundle
        """
        return self.homogeneous_space

    def chern_classes(self):
        r"""
        Return the list of homogeneous parts of Chern classes of this vector bundle
        """

        def class_from_weight(weight):
            return sum(
                weight[i] * self.homogeneous_space.x[i]
                for i in range(
                    self.homogeneous_space.parabolic_subgroup.ambient_space_dimension()
                )
            )

        cc = prod(
            (1 + class_from_weight(vector(w))) ^ i
            for w, i in self.weight_multiplicities.items()
        )

        return [homogeneous_part(cc, i) for i in ((0.0).self.homogeneous_space.dim)]


class IrreducibleEquivariantVectorBundle(EquivariantVectorBundle):
    def __init__(self, homogeneous_space, weight) -> None:
        r"""

        Constructor of this class

        This constructor takes the base homogeneous space `G/P` and weights of `G` with their multiplicities. This construct an equivariant vector bundle associated to the highest weight representation of `G` associated to the weight ``weight``.

        INPUT:

        - ``weight`` -- list of integers -- this represents the coefficients of the fundamental weights

        """
        self.weight = weight
        super().__init__(
            homogeneous_space,
            homogeneous_space.parabolic_subgroup.weight_multiplicities(weight),
        )

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.homogeneous_space} associated to {self.weight}"
