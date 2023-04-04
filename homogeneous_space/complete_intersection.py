from sage.all import *
from homogeneous_space.homogeneous_space import (
    EquivariantVectorBundle,
    HomogeneousSpace,
)
from homogeneous_space.interfaces import IVariety

homogeneous_part = lambda F, degree: sum(
    c * m for c, m in F if m.total_degree() == degree
)


class CompleteIntersection(IVariety):
    def __init__(
        self, flag_variety: HomogeneousSpace, vector_bundle: EquivariantVectorBundle
    ) -> None:
        self.flag_variety = flag_variety
        self.vector_bundle = vector_bundle
        self.dim = self.flag_variety.dimension() - self.vector_bundle.rank()

    def __repr__(self) -> str:
        return (
            f"a complete intersection of {self.flag_variety} and {self.vector_bundle}"
        )

    def dimension(self) -> int:
        return self.dim

    def tangent_bundle(self):
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
        def class_from_weight(weight):
            return sum(
                weight[i] * self.flag_variety.x[i]
                for i in range(
                    self.flag_variety.parabolic_subgroup.ambient_space_dimension()
                )
            )

        def geometric_sequence(n, x):
            return sum(x ^ i for i in ((0.0).n))

        cc = prod(1 + x for x in self.flag_variety.tangent_weights) * prod(
            geometric_sequence(self.dim, -class_from_weight(vector(w))) ^ i
            for w, i in self.vector_bundle.weight_muliplicities.items()
        )

        return [homogeneous_part(cc, i) for i in ((0.0).self.dim)]

    def numerical_integration_by_localization(self, f):
        top_of_f = homogeneous_part(f, self.dim)
        c_top = (
            self.vector_bundle.chern_classes()[self.vector_bundle.rank()]
            if self.dim >= 0
            else 0
        )
        return self.flag_variety.numerical_integration_by_localization(top_of_f * c_top)

    def integration(self, f) -> int:
        return self.numerical_integration_by_localization(f)
