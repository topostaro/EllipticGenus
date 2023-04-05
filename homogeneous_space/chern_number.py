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

from sage.all import prod
from homogeneous_space.interfaces import IVariety

# `degree`次部分を取り出す関数
homogeneous_part = lambda F, degree: sum(
    c * m for c, m in F if m.total_degree() == degree
)


def chern_number(variety: IVariety, degrees: list) -> int:
    if sum(d for d in degrees) != variety.dimension():
        return 0
    else:
        chern_classes = variety.chern_classes()
        return variety.integration(prod([chern_classes[d] for d in degrees]))
