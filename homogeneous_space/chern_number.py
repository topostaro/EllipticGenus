r"""
Computation of Chern numbers of varieties
================================================

This module implements a computation of Chern numbers of varieties.

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
    r"""

    Return the Chern number of ``variety`` with degree ``degree``

    INPUT:

    - ``variety`` -- ``IVariety``

    - ``degree`` -- list of integers -- this list represents the integrant. For example, ``[1,1,3]`` represents `c_1 c_1 c_3`.

    OUTPUT:

    the Chern number of ``variety`` with degree ``degree``

    """
    if sum(d for d in degrees) != variety.dimension():
        return 0
    else:
        chern_classes = variety.chern_classes()
        return variety.integration(prod([chern_classes[d] for d in degrees]))
