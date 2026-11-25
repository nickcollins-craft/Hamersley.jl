#=
Kinematic nodes.

These nodes build the standard kinematic quantities of continuum mechanics from a
displacement field. They are nodes in their own right (rather than constitutive inputs
computed somewhere else) because their linearisation is where the geometric stiffness of
a finite-strain problem comes from.
=#


"""
    AbstractKinematicNode

Nodes that compute a kinematic quantity from a displacement field.
"""
abstract type AbstractKinematicNode <: AbstractNode end


"""
    DeformationGradient(field)

The deformation gradient ``F = I + ∇_0 u`` of the displacement field `field`.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct DeformationGradient <: AbstractKinematicNode
    field::Symbol
end


"""
    GreenLagrangeStrain(field)

The Green–Lagrange strain ``E = \\tfrac{1}{2}(F^T F - I)`` of the displacement field
`field`.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct GreenLagrangeStrain <: AbstractKinematicNode
    field::Symbol
end


"""
    InfinitesimalStrain(field)

The infinitesimal strain ``ε = \\operatorname{sym}(∇u)`` of the displacement field
`field`. It is only meaningful under [`SmallStrain`](@ref).

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct InfinitesimalStrain <: AbstractKinematicNode
    field::Symbol
end


"""
    VolumeRatio(field)

The volume ratio ``J = \\det F`` of the displacement field `field`. Under small strains it
is linearised to ``J ≈ 1 + \\operatorname{tr} ε``.

# Fields
- `field::Symbol`: the name of the displacement field ``u``.
"""
struct VolumeRatio <: AbstractKinematicNode
    field::Symbol
end


operand_field(node::AbstractKinematicNode) = node.field

# The quantity tag of each kinematic node does not depend on the field's declaration, but
# we still check that the field exists so that undeclared fields give no tag.
quantity(node::DeformationGradient, fields::Tuple) =
    field_quantity(_ -> DeformationGradientTensor(), node, fields)
quantity(node::GreenLagrangeStrain, fields::Tuple) =
    field_quantity(_ -> GreenLagrangeStrainTensor(), node, fields)
quantity(node::InfinitesimalStrain, fields::Tuple) =
    field_quantity(_ -> InfinitesimalStrainTensor(), node, fields)
quantity(node::VolumeRatio, fields::Tuple) =
    field_quantity(_ -> JacobianDeterminant(), node, fields)

node_fields(node::AbstractKinematicNode, ::AbstractKinematicAssumption) =
    [operand_field(node)]

subnodes(::AbstractKinematicNode) = ()


"""
    is_affine(node, assumption)

Return `true` when the value of `node` is an affine function of the fields it depends on
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
