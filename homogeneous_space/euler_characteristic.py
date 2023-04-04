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
