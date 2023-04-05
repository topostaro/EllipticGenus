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

from homogeneous_space.interfaces import IVariety, IVectorBundle


def euler_characteristic(variety: IVariety, vector_bundle: IVectorBundle) -> int:
    chern_character = vector_bundle.chern_character()
    todd_classes = variety.todd_classes()

    return variety.integration(
        sum(
            chern_character[i] * todd_classes[variety.dimension() - i]
            for i in ((0.0).variety.dimension())
        )
    )
