# EllipticGenus

A Sagemath package for computing elliptic genera of homogeneous spaces and complete intersections in them, and representing elliptice genera of any complex manifold by using Chern numbers.

##  Structure of this module

- ellptic_genus
    This module implements the core functions of computing elliptic genera. 
    - elliptic_genus.py
        - `elliptic_genus(variety, k)` computes the elliptic genera of `variety` up to degree `k`. The argument `variety` should be an object of `IVariety`, equivalently, the object must contain the data of Chern numbers.
        - `elliptic_genus_chernnum(dim, k)` outputs the elliptic genus of varieties of dimension `dim` up to degree `k` represented symbolically by Chern numbers.
    - utils.py
        This file contains private functions.

- weak_Jacobi_form
    This module implements a computation of a basis of weak Jacobi forms of weight `0` and indices that are either integral or half-integral. 
    - basis.py
