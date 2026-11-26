#=
Kinematic vertices.

These vertices build the standard kinematic quantities of continuum mechanics from a
displacement field. They are vertices in their own right (rather than constitutive inputs
computed somewhere else) because their linearisation is where the geometric stiffness of
a finite-strain problem comes from.
=#


"""
    AbstractKinematicVertex

Vertices that compute a kinematic quantity from a displacement field.
"""
abstract type AbstractKinematicVertex <: AbstractVertex end


"""
    DeformationGradient(field)

The deformation gradient ``F = I + ∇_0 u`` of the displacement field `field`.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct DeformationGradient <: AbstractKinematicVertex
    field::Symbol
end


"""
    GreenLagrangeStrain(field)

The Green–Lagrange strain ``E = \\tfrac{1}{2}(F^T F - I)`` of the displacement field
`field`.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct GreenLagrangeStrain <: AbstractKinematicVertex
    field::Symbol
end


"""
    InfinitesimalStrain(field)

The infinitesimal strain ``ε = \\operatorname{sym}(∇u)`` of the displacement field
`field`. It is only meaningful under [`SmallStrain`](@ref).

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct InfinitesimalStrain <: AbstractKinematicVertex
    field::Symbol
end


"""
    VolumeRatio(field)

The volume ratio ``J = \\det F`` of the displacement field `field`. Under small strains it
is linearised to ``J ≈ 1 + \\operatorname{tr} ε``.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct VolumeRatio <: AbstractKinematicVertex
    field::Symbol
end


operand_field(vertex::AbstractKinematicVertex) = vertex.field

# The quantity tag of each kinematic vertex does not depend on the field's declaration, but
# we still check that the field exists so that undeclared fields give no tag.
quantity(vertex::DeformationGradient, fields::Tuple) =
    field_quantity(_ -> DeformationGradientTensor(), vertex, fields)
quantity(vertex::GreenLagrangeStrain, fields::Tuple) =
    field_quantity(_ -> GreenLagrangeStrainTensor(), vertex, fields)
quantity(vertex::InfinitesimalStrain, fields::Tuple) =
    field_quantity(_ -> InfinitesimalStrainTensor(), vertex, fields)
quantity(vertex::VolumeRatio, fields::Tuple) =
    field_quantity(_ -> JacobianDeterminant(), vertex, fields)

vertex_fields(vertex::AbstractKinematicVertex, ::AbstractKinematicAssumption) =
    [operand_field(vertex)]

subvertices(::AbstractKinematicVertex) = ()


"""
    is_affine(vertex, assumption)

Return `true` when the value of `vertex` is an affine function of the fields it depends on
under `assumption`.

``F = I + ∇_0 u`` and ``ε`` are affine in ``u``. ``E`` and ``J`` are not under finite
strains, but under small strains they are linearised (to ``ε`` and ``1 + \\operatorname{tr}
ε``) and so are affine.
"""
is_affine(::DeformationGradient, ::AbstractKinematicAssumption) = true
is_affine(::InfinitesimalStrain, ::AbstractKinematicAssumption) = true
is_affine(::GreenLagrangeStrain, ::SmallStrain) = true
is_affine(::GreenLagrangeStrain, ::FiniteStrain) = false
is_affine(::VolumeRatio, ::SmallStrain) = true
is_affine(::VolumeRatio, ::FiniteStrain) = false
