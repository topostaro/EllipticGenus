from sage.rings.rational_field import Q
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.lazy_series_ring import LazyLaurentSeriesRing


R0 = LaurentPolynomialRing(Q, "y")
R = LazyLaurentSeriesRing(R0, "q")
