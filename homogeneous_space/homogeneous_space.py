from sage.all import *
from parabolic import ParabolicSubgroup
from interfaces import *


class HomogeneousSpace(IVariety):
    def __init__(self, parabolic_subgroup) -> None:
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
        return f"a flag_variety associated to {self.parabolic_subgroup}"

    def dimension(self) -> int:
        return self.dim

    def tangent_bundle(self) -> IVectorBundle:
        tangent_weights = {
            w: 1
            for w in set(self.parabolic_subgroup.R_G.positive_roots())
            - set(self.parabolic_subgroup.positive_roots())
        }
        return EquivariantVectorBundle(self, tangent_weights)

    # 次数ごとのchern類
    def chern_classes(self):
        return [
            homogeneous_part(prod(1 + x for x in self.tangent_weights), i)
            for i in ((0.0).self.dim)
        ]

    def numerical_integration_by_localization(self, f):
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
        return self.numerical_integration_by_localization(f)


class EquivariantVectorBundle(IVectorBundle):
    def __init__(self, flag_variety, weight_muliplicities) -> None:
        self.flag_variety = flag_variety
        self.weight_muliplicities = weight_muliplicities
        self.rk = sum(v for v in self.weight_muliplicities.values())

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.flag_variety} associated to {self.weight_muliplicities}"

    def rank(self) -> int:
        return self.rk

    def base(self) -> HomogeneousSpace:
        return self.flag_variety

    def chern_classes(self):
        def class_from_weight(weight):
            return sum(
                weight[i] * self.flag_variety.x[i]
                for i in range(
                    self.flag_variety.parabolic_subgroup.ambient_space_dimension()
                )
            )

        cc = prod(
            (1 + class_from_weight(vector(w))) ^ i
            for w, i in self.weight_muliplicities.items()
        )

        return [homogeneous_part(cc, i) for i in ((0.0).self.flag_variety.dim)]


class IrreducibleEquivariantVectorBundle(EquivariantVectorBundle):
    def __init__(self, flag_variety, weight) -> None:
        self.weight = weight
        super().__init__(
            flag_variety, flag_variety.parabolic_subgroup.weight_muliplicities(weight)
        )

    def __repr__(self) -> str:
        return f"an equivariant vector bundle on {self.flag_variety} associated to {self.weight}"
