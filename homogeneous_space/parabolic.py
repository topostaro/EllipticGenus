r"""
Class representing a parabolic subgroup representing a parabolic subgroup of a reductive algebraic group
================================================

This module contains a class representing a parabolic subgroup representing a parabolic subgroup of a reductive algebraic group.
The class provides methods to compute its simple roots, positive roots, and the weight set of its highest weight representations.

EXAMPLE:

    sage: from homogeneous_space.parabolic import ParabolicSubgroup
    sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
    sage: P
    the parabolic subgroup of ['A', 3] with crossed-out nodes [1]
    sage: P.dynkin_diagram()
    X---O---O
    1   2   3
    A3 with node 1 marked
    sage: P.positive_roots()
    [(0, 1, -1, 0), (0, 0, 1, -1), (0, 1, 0, -1)]
    sage: P.weight_multiplicities((1, 2, 0))
    {(3, 1, 1, 0): 1,
     (3, 1, 0, 1): 1,
     (3, 2, 0, 0): 1,
     (3, 0, 1, 1): 1,
     (3, 0, 2, 0): 1,
     (3, 0, 0, 2): 1}

AUTHORS:

- KENTA KOBAYASHI and AKIHITO NAKAMURA (2023-04-04): initial version

REFERENCES:

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

from sage.all import matrix, vector, WeylCharacterRing, CartanType


# 最高ウェイト表現のウェイトを最高ウェイトとの差で表す関数。
def root_difference_multiplicities(character_ring, highest_weight) -> dict:
    r"""

    Return the dictionary which represents the weights of the highest weight representation in the character ring.

    INPUT:

    - ``character_ring`` -- WeylCharacterRing

    - ``highest_weight`` -- weight of the ``character_ring``

    OUTPUT:

    the dictionary, where the keys are tuples of integers and the values are the multiplicities of the weight which is the sum of ``highest_weight`` and simple roots with the coefficients of their keys.


    EXAMPLE:

        sage: from homogeneous_space.parabolic import root_difference_multiplicities
        sage: root_difference_multiplicities(WeylCharacterRing("A1"), (2,1))
        {(-1,): 1, (0,): 1}

    """
    weight_multiplicities = character_ring(highest_weight).weight_multiplicities()
    A = matrix([vector(sr) for sr in character_ring.simple_roots()]).transpose()

    result = {}

    for k, v in weight_multiplicities.items():
        Y = vector(k) - vector(highest_weight)
        result[tuple(A.solve_right(Y))] = v

    return result


# Parabolic subgroupを表すクラス
class ParabolicSubgroup:
    r"""

    Class representing a parabolic subgroup of a reductive algebraic group

    This class contains the data of crossed Dynkin diagrams. This provides methods to compute its simple roots, positive roots, and the weight set of its highest weight representations.

    """

    def __init__(self, G, L, crossed_out_nodes) -> None:
        r"""

        Constructor of this class

        INPUT:

        - ``G`` -- CartanType -- the Cartan type of the parent Lie group

        - ``L`` -- CartanType | None -- the Cartan type of the Levi subgroup or None if the CartanType is empty

        - crossed_out_nodes -- list of integers -- the indices of crossed-out nodes


        OUTPUT:

        the initialized object

        EXAMPLE:

            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
            sage: P
            the parabolic subgroup of ['A', 3] with crossed-out nodes [1]

        """
        self.G = G
        self.R_G = WeylCharacterRing(self.G)

        self.L = L

        if L != None:
            self.R_L = WeylCharacterRing(self.L)
        else:
            self.R_L = None

        self.crossed_out_nodes = crossed_out_nodes

    def __repr__(self) -> str:
        return f"the parabolic subgroup of {self.G} with crossed-out nodes {self.crossed_out_nodes}"

    def dynkin_diagram(self):
        r"""

        Return the crossed Dynkin diagram of the parabolic subgroup

        EXAMPLE:

            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
            sage: P.dynkin_diagram()
            X---O---O
            1   2   3
            A3 with node 1 marked

        """
        return self.G.marked_nodes(self.crossed_out_nodes).dynkin_diagram()

    def ambient_space(self):
        r"""

        Return the ambient space of the parent group of the parabolic subgroup

        """
        return self.R_G.space()

    # Sagemathのライブラリでは1-indexを用いるが、ここでは簡単のためとりあえず0-indexで表す
    def simple_roots(self):
        r"""

        Return the simple root of the parabolic subgroup

        EXAMPLE:

            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
            sage: P.simple_roots()
            [(0, 1, -1, 0), (0, 0, 1, -1)]

        """

        return [
            self.R_G.simple_roots()[i]
            for i in set(range(1, self.G.rank() + 1)) - set(self.crossed_out_nodes)
        ]

    def positive_roots(self):
        r"""

        Return the positive root of the parabolic subgroup

        EXAMPLE:

            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
            sage: P.positive_roots()
            [(0, 1, -1, 0), (0, 0, 1, -1), (0, 1, 0, -1)]

        """

        L = self.ambient_space()

        # positive rootからcutoutされるsimple rootを引き, positiveでなければuncrossed nodeから生成されると判定
        roots = [pr for pr in self.G.root_system().root_lattice().positive_roots()]
        for i in self.crossed_out_nodes:
            roots = [
                pr
                for pr in roots
                if not (
                    pr - (self.G.root_system().root_lattice().simple_roots())[i]
                ).is_positive_root()
            ]

        # ambient spaceの元に変換
        return [L(pr) for pr in roots]

    # 引数の`weight`は基本ウェイトを基底にして表示したもの
    def weight_multiplicities(self, weight) -> dict:
        r"""

        Return the dictionary, where the keys are weights and the values are their multiplicities

        INPUT:

        - ``weight`` -- list of integers -- the list of coefficients of fundamental weights

        EXAMPLE:

            sage: from homogeneous_space.parabolic import ParabolicSubgroup
            sage: P = ParabolicSubgroup(CartanType('A3'), CartanType('A2'), [1])
            sage: P.weight_multiplicities((1, 2, 0))
            {(3, 1, 1, 0): 1,
             (3, 1, 0, 1): 1,
             (3, 2, 0, 0): 1,
             (3, 0, 1, 1): 1,
             (3, 0, 2, 0): 1,
             (3, 0, 0, 2): 1}

        """
        if self.L == None:
            fws_G = [fw for fw in self.R_G.fundamental_weights()]
            weight_for_G = sum(weight[i] * fws_G[i] for i in range(self.G.rank()))

            return {weight_for_G: 1}

        # GとLのディンキン図の頂点のずれを補正する関数
        def correct_index(index: int) -> int:
            for i in range(len(self.crossed_out_nodes)):
                if index + i < self.crossed_out_nodes[i]:
                    return index + i
            return index + len(self.crossed_out_nodes)

        fws_L = [
            fw for fw in self.R_L.fundamental_weights()
        ]  # conversion from 1-index to 0-index
        weight_for_L = [
            weight[i - 1]
            for i in set(range(1, len(weight) + 1)) - set(self.crossed_out_nodes)
        ]
        weight_for_L = sum(weight_for_L[i] * fws_L[i] for i in range(self.L.rank()))

        fws_G = [fw for fw in self.R_G.fundamental_weights()]
        weight_for_G = sum(weight[i] * fws_G[i] for i in range(self.G.rank()))

        mul_set = root_difference_multiplicities(self.R_L, weight_for_L)

        result = {}
        for k, v in mul_set.items():
            w = weight_for_G + sum(
                k[i - 1] * self.R_G.simple_roots()[correct_index(i)]
                for i in range(1, self.L.rank() + 1)
            )
            result[w] = v

        return result
