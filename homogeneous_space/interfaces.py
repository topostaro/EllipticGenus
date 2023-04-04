r"""
Abstract classes of variety and vector bundles and some operations for them
================================================

This module contains abstract classes:
    - ``IVariety`` -- an interface of varieties,
    - ``IVectorBundle`` -- an interface of vector bundles.
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


from sage.all import *
from abc import ABC, abstractmethod
from parabolic import ParabolicSubgroup


class IVectorBundle(ABC):
    r"""

    Dummy definition for type annotations

    """
    pass


class IVariety(ABC):
    r"""

    Abstract class of varieties

    """

    @abstractmethod
    def dimension(self) -> int:
        r"""
        Return the dimension of this variety
        """
        pass

    @abstractmethod
    def tangent_bundle(self) -> IVectorBundle:
        r"""
        Return the tangent bundle of this variety
        """
        pass

    def cotangent_bundle(self) -> IVectorBundle:
        r"""
        Return the cotangent bundle of this variety
        """
        return self.tangent_bundle().dual()

    @abstractmethod
    def integration(self, f) -> int:
        r"""
        Return the integration value of the cohomology class ``f`` on this variety

        INPUT:

        - ``f`` -- a cohomology class

        OUTPUT:

        the integration of ``f`` on this variety

        """
        pass

    def chern_classes(self) -> list:
        r"""
        Return the degreewise list of Chern classes of the tangent bundle of this variety
        """
        return self.tangent_bundle().chern_classes()

    def todd_classes(self) -> list:
        r"""
        Return the todd classes of the tangent bundle of this variety
        """
        return self.tangent_bundle().todd_classes()


import re

singular.lib("chern.lib")


class IVectorBundle(ABC):
    r"""

    Abstract class of vector bundle

    """

    @abstractmethod
    def base(self) -> IVariety:
        r"""
        Return the base space of this vector bundle
        """
        pass

    @abstractmethod
    def rank(self) -> int:
        r"""
        Return the rank of this vector bundle
        """
        pass

    @abstractmethod
    def chern_classes(self) -> list:
        r"""
        Return the degreewise list of Chern classes of this vector bundle
        """
        pass

    # 演算子のオーバーロード
    def __add__(self, other) -> IVectorBundle:
        r"""

        Return the direct sum of vector bundles

        """
        return direct_sum(self, other)

    def __mul__(self, other) -> IVectorBundle:
        r"""

        Return the tensor product of vector bundles

        """
        return tensor_product(self, other)

    def chern_character(self) -> list:
        r"""
        Return the list of homogeneous parts of Chern characters of this vector bundle
        """
        len_cc = len(self.chern_classes()) - 1

        ring_for_ch = PolynomialRing(
            QQ,
            [f"c{i}_E" for i in ((1.0).len_cc)],
            order=TermOrder("wdeglex", tuple((1.0).len_cc)),
        )

        # Using Singular, calclate universal formula of chern character
        singular.lib("chern.lib")
        r = singular.ring(0, f"(c(1..{len_cc}))", "dp")
        l = singular.list(f"c(1..{len_cc})")
        ch_str_list = singular.chAll(
            l, self.base().dimension()
        ).sage_structured_str_list()
        chern_character = [ring_for_ch(self.rank())] + [
            ring_for_ch(re.sub(r"c\(([0-9]+)\)", r"c\1_E", s)) for s in ch_str_list
        ]

        chern_classes = self.chern_classes()[1:]

        return [
            chern_character[i](chern_classes) for i in ((0.0).self.base().dimension())
        ]

    def chern_character_total(self):
        r"""
        Return the Chern characters of this vector bundle
        """
        return sum(c for c in self.chern_character())

    def todd_classes(self) -> list:
        r"""
        Return the list of homogeneous parts of Todd classes of this vector bundle
        """
        len_cc = len(self.chern_classes()) - 1

        ring_for_td = PolynomialRing(
            QQ,
            [f"c{i}_M" for i in ((1.0).len_cc)],
            order=TermOrder("wdeglex", tuple((1.0).len_cc)),
        )

        # Using Singular, calculate universal formula of Todd classes
        r = singular.ring(0, f"(c(1..{len_cc}))", "dp")
        l = singular.list(f"c(1..{len_cc})")
        todd_str_list = singular.todd(l).sage_structured_str_list()
        todd_classes = [ring_for_td(1)] + [
            ring_for_td(re.sub(r"c\(([0-9]+)\)", r"c\1_M", s)) for s in todd_str_list
        ]

        chern_classes = self.chern_classes()[1:]

        return [todd_classes[i](chern_classes) for i in ((0.0).self.base().dimension())]

    def dual(self) -> IVectorBundle:
        r"""
        Return the dual vector bundle of this vector bundle
        """
        vector_bundle = self

        class VB(IVectorBundle):
            def rank(self) -> int:
                return vector_bundle.rank()

            def base(self) -> IVariety:
                return vector_bundle.base()

            def chern_classes(self) -> list:
                chern_classes = vector_bundle.chern_classes()
                return [(-1) ^ i * chern_classes[i] for i in range(len(chern_classes))]

            def __repr__(self) -> str:
                return f"the dual vector bundle of {vector_bundle}"

        return VB()


def direct_sum(vector_bundle1: IVectorBundle, vector_bundle2: IVectorBundle):
    r"""

    Return the direct sum of vector bundles

    """
    if vector_bundle1.base() != vector_bundle2.base():
        raise TypeError("Not match bases of vector bundles")

    rank = vector_bundle1.rank() + vector_bundle2.rank()
    base = vector_bundle1.base()

    cc = sum(c1 for c1 in vector_bundle1.chern_classes()) * sum(
        c2 for c2 in vector_bundle2.chern_classes()
    )
    chern_classes = [
        homogeneous_part(cc, i) for i in ((0.0).vector_bundle1.base().dimension())
    ]

    class VB(IVectorBundle):
        def rank(self) -> int:
            return rank

        def base(self) -> IVariety:
            return base

        def chern_classes(self) -> list:
            return chern_classes

        def __repr__(self) -> str:
            return f"the direct sum of {vector_bundle1} and {vector_bundle2}"

    return VB()


def tensor_product(vector_bundle1: IVectorBundle, vector_bundle2: IVectorBundle):
    r"""

    Return the tensor product of vector bundles

    """
    if vector_bundle1.base() != vector_bundle2.base():
        raise TypeError("Not match bases of vector bundles")

    rank1 = vector_bundle1.rank()
    rank2 = vector_bundle2.rank()

    len_cc1 = len(vector_bundle1.chern_classes()) - 1
    len_cc2 = len(vector_bundle2.chern_classes()) - 1

    ring_for_Es = PolynomialRing(
        QQ,
        [f"c{i}_E1" for i in ((1.0).len_cc1)] + [f"c{i}_E2" for i in ((1.0).len_cc2)],
    )

    r = singular.ring(0, f"(c(1..{len_cc1}), C(1..{len_cc2}))", "dp")
    l1 = singular.list(f"c(1..{len_cc1})")
    l2 = singular.list(f"C(1..{len_cc2})")

    ch_str_list = singular.chProd(rank1, l1, rank2, l2).sage_structured_str_list()
    ch_prod = [
        re.sub(r"c\(([0-9]+)\)", r"c\1_E1", s)
        for s in ch_str_list[: vector_bundle1.base().dimension()]
    ]
    ch_prod = [ring_for_Es(re.sub(r"C\(([0-9]+)\)", r"c\1_E2", s)) for s in ch_prod]

    cc = [1] + [
        cp(vector_bundle1.chern_classes()[1:] + vector_bundle2.chern_classes()[1:])
        for cp in ch_prod
    ]

    rank = vector_bundle1.rank() * vector_bundle2.rank()
    base = vector_bundle1.base()

    class VB(IVectorBundle):
        def rank(self) -> int:
            return rank

        def base(self) -> IVariety:
            return base

        def chern_classes(self) -> list:
            return cc

        def __repr__(self) -> str:
            return f"the tensor product of {vector_bundle1} and {vector_bundle2}"

    return VB()
