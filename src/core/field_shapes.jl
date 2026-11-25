#=
Field shapes.

The shape of a field says what kind of value it takes at a point (a scalar, a vector, or a
tensor). Shapes are needed both by fields themselves and by quantity tags such as the
gradient of a scalar versus the gradient of a vector, which have different conjugates.
=#


"""
    AbstractFieldShape

Supertype of the value shapes a field (or a quantity derived from one) can take.
"""
abstract type AbstractFieldShape end


"""
    ScalarField()

A scalar-valued field, such as the temperature ``θ`` or the pore pressure ``p``.
"""
struct ScalarField <: AbstractFieldShape end


"""
    VectorField()

A vector-valued field, such as the displacement ``u`` or a flux ``q``. Its dimension is the
spatial dimension of the mesh, which is only known once the grid is.
"""
struct VectorField <: AbstractFieldShape end


"""
    SymmetricTensorField()

A symmetric second-order tensor field, such as an interpolated Cauchy stress ``σ`` in a
stress-based mixed formulation.
"""
struct SymmetricTensorField <: AbstractFieldShape end


"""
    TensorField()

A general second-order tensor field, such as the relaxed displacement gradient ``ψ ≈ ∇u``
of a micromorphic formulation.
"""
struct TensorField <: AbstractFieldShape end
