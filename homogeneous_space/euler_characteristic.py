r"""
Computation of Euler characteristic of complex vector bundles
================================================

This module implements a computation of Euler characteristic of complex vector bundles.

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
    r"""

    Return the Euler characteristic of ``vector_bundle`` on ``variety``

    INPUT:

    - ``variety`` -- object of ``IVariety`` -- the base space.

    - ``vector_bundle`` -- object of ``IVectorBundle``

    OUTPUT:

    the Euler characteristic of ``vector_bundle`` on ``variety``, that is computed as the integration of the multiplication of the Chern character of ``vector_bundle`` and the Todd genus of ``variety``.

    ..MATH::

    \chi (X, E) = \int_X ch(E) td(X)

    """
    chern_character = vector_bundle.chern_character()
    todd_classes = variety.todd_classes()

    return variety.integration(
        sum(
            chern_character[i] * todd_classes[variety.dimension() - i]
            for i in range(0, variety.dimension() + 1)
        )
    )
