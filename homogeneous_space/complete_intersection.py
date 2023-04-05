r"""
Class of varieties which is the zeros of a general section of an equivariant vector bundle of a homogeneous space
================================================

This module contains a class:
    - ``CompleteIntersection`` -- an implementation of ``IVariety``, representing the zeros of a general section of an equivariant vector bundle of a homogeneous space
This is specialized in computing Chern characters and Todd classes from Chern classes.

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

from sage.all import prod, vector
from homogeneous_space.homogeneous_space import (
    EquivariantVectorBundle,
    HomogeneousSpace,
)
from homogeneous_space.interfaces import IVariety, IVectorBundle

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

        """
        self.homogeneous_space = homogeneous_space
        self.vector_bundle = vector_bundle
        self.dim = self.homogeneous_space.dimension() - self.vector_bundle.rank()

    def __repr__(self) -> str:
        return f"a complete intersection of {self.homogeneous_space} and {self.vector_bundle}"

    def dimension(self) -> int:
        r"""
        Return the dimension of this variety
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
        """

        def class_from_weight(weight):
            return sum(
                weight[i] * self.homogeneous_space.x[i]
                for i in range(
                    self.homogeneous_space.parabolic_subgroup.ambient_space_dimension()
                )
            )

        def geometric_sequence(n, x):
            return sum(x ^ i for i in ((0.0).n))

        cc = prod(1 + x for x in self.homogeneous_space.tangent_weights) * prod(
            geometric_sequence(self.dim, -class_from_weight(vector(w))) ^ i
            for w, i in self.vector_bundle.weight_multiplicities.items()
        )

        return [homogeneous_part(cc, i) for i in ((0.0).self.dim)]

    def numerical_integration_by_localization(self, f):
        r"""

        Return the numerical computation of the integration of equivariant cohomology classes.

        INPUT:

        - ``f`` -- an equivariant cohomology class on this variety

        OUTPUT:

        the numerical computation of the integration of the  equivariant cohomology class ``f``.

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
