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
