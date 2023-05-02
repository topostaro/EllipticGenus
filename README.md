# EllipticGenus

A SageMath package for computing elliptic genera of homogeneous spaces and complete intersections in them, and representing elliptic genera of any complex manifold by using Chern numbers.

##  Structure of this module

- elliptic_genus

    This module implements the core functions of computing elliptic genera of almost complex manifolds.

    - elliptic_genus.py

        - `elliptic_genus(manifold, k)` computes the elliptic genera of `manifold` with the terms of `q`-variable up to degree `k`. The argument `manifold` should be an object of `AlmostComplexManifold`, equivalently, the object must contain the data of Chern numbers.

        - `elliptic_genus_chernnum(dim, k)` computes the elliptic genus of a almost complex manifold of dimension `dim` with the terms of `q`-variable up to degree `k` represented symbolically by Chern numbers.

    - utils.py

        This file contains private functions.

- weak_Jacobi_form

    This module implements a computation of a basis of weak Jacobi forms of weight 0 and indices that are either integral or half-integral. 

    - basis.py

        - `basis_integral(index)` returns a basis of the space of weak Jacobi forms of weight 0 and index `index`. `index` should be a non-negative integer.

        - `basis_half_integral(double_index)` returns a basis of the space of weak Jacobi forms of weight 0 and index `double_index / 2` for even `double_index`, otherwise return a list of a basis of the space of weak Jacobi forms of weight 0 and index `double_index / 2` multiplied by $y^{1/2}$.

    - eisenstein.py

        - `eisenstein(k)` computes the normalized Eisenstein series of weight `k`.

- homogeneous_space

    This module provides various functions which compute the integrations of equivariant cohomology classes on homogeneous spaces and complete intersections of them. For example, there are functions for computation of Euler characteristic of equivariant vector bundles and for computation of Chern numbers of such manifolds.

    - parabolic.py

        - The class `ParabolicSubgroup` represents a parabolic subgroup of a reductive algebraic group. This provides methods to compute its simple roots, positive roots, and the weight set of its highest weight representations.

    - interfaces.py

        - The abstract class `AlmostComplexManifold` represents an almost complex manifold which has the methods providing Chern classes and their integrations.

        - The abstract class `VectorBundle` represents a vector bundle which has the methods providing Chern classes. The operators `+` and `*` are overloaded by the following functions `direct_sum` and `tensor_product` respectively.

        - `direct_sum(vector_bundle1, vector_bundle2)` returns the direct sum of vector bundles as an object implementing `VectorBundle`.

        - `tensor_product(vector_bundle1, vector_bundle2)` returns the direct sum of vector bundles as an object implementing `VectorBundle`.

    - homogeneous_space.py

        - The class `HomogeneousSpace` represents a homogeneous space, which is the quotient by a parabolic subgroup. The integration of cohomology classes is implemented by using numerical computation of localization.

        - The class `EquivariantVectorBundle` represents a equivariant vector bundle on a homogeneous space. This provides the Chern classes of the vector bundle.

        - The class `IrreducibleEquivariantVectorBundle` is a variation of `EquivariantVectorBundle`, which is specialized for irreducible ones.

    -  complete_intersection.py

        - The class `CompleteIntersection` represents the zeros of a general section of an equivariant vector bundle of a homogeneous space.

    - euler_characteristic.py

        - `euler_characteristic(manifold, vector_bundle)` computes the Euler characteristic of `vector_bundle` on `manifold` by using the method `integration()` of `manifold`.

    - chern_number.py

        - `chern_number(manifold, degrees)` provides the Chern number of `manifold` with degree `degree`. The argument `degree` should be a list of positive integers. For example, `[1,1,3]` represents $c_1 c_1 c_3$.

## Installation

This package can be installed using pip in SageMath.

```
sage --pip install git+https://github.com/topostaro/EllipticGenus.git
```

## Usage

We give an example to show how to use this package. We compute the elliptic genus of quintic Calabi--Yau 3-fold, and more.

1. First, construct quintic 3-fold.

```
sage: from homogeneous_space import *
sage: P = ParabolicSubgroup(CartanType('A4'), CartanType('A3'), [1])
sage: Proj4 = HomogeneousSpace(P)
sage: L = IrreducibleEquivariantVectorBundle(Proj4, (5, 0, 0, 0, 0))
sage: Quintic = CompleteIntersection(L)
```
2. Compute the elliptic genus.

```
sage: from elliptic_genus import elliptic_genus
sage: elliptic_genus(Quintic, 2)
-100*y - 100*y^2 + (100*y^-1 - 100*y - 100*y^2 + 100*y^4)*q + (100*y^-2 + 100*y^-1 - 200*y - 200*y^2 + 100*y^4 + 100*y^5)*q^2 + O(q^3)
```

3. Compute the Chern number $\int c_3$ of the quintic 3-fold.

```
sage: chern_number(Quintic, [3])
-200
```
4. Compute the Euler characteristic of the tangent bundle of the quintic 3-fold.

```
sage: euler_characteristic(Quintic, Quintic.tangent_bundle())
-100
```