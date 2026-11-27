#=
Boundary vertices.

Natural and Robin conditions are compiled into boundary terms of the weak form, so that
the linearity, coupling and measure analyses (and, later, assembly) treat them exactly as
they treat every other term. These are the vertices those terms are built from: the trace
of the test function on the boundary, a prescribed load, and the flux of a Robin
condition.

A dead load is constant, so its term is affine and depends on no field. A follower load is
integrated over the current boundary, whose pull-back (Nanson's formula) depends on the
motion; that dependence belongs to the term's measure, not to the load. A load whose value
depends on the current position (a hydrostatic pressure on a deforming body) depends on the
motion through its value as well.
=#


"""
    BoundaryTrace(field, trace)

The trace `trace` of `field` on the boundary: the test side of the boundary term of a
natural condition whose trace is not the value (``∂_{n_0} δu`` for a double force,
``δq \\cdot n_0`` for a pressure on a flux field).

# Fields
- `field::Symbol`: the name of the field.
- `trace::T`: the [`AbstractTrace`](@ref).
"""
struct BoundaryTrace{T<:AbstractTrace} <: AbstractFieldOperator
    field::Symbol
    trace::T
end


# Only the value trace carries a tag of its own: the field's.
quantity(operator::BoundaryTrace{ValueTrace}, fields::Tuple) =
    field_quantity(field_value_quantity, operator, fields)
quantity(::BoundaryTrace, ::Tuple) = nothing


"""
    BoundaryLoad(quantity, value, components)

The load of a natural condition: the boundary quantity `quantity`, with the prescribed
value `value`, on the components `components`.

# Fields
- `quantity::Q`: the [`AbstractBoundaryQuantity`](@ref) tag.
- `value::P`: the [`AbstractPrescription`](@ref).
- `components::C`: the [`AbstractComponentSelection`](@ref).
"""
struct BoundaryLoad{
    Q<:AbstractBoundaryQuantity,
    P<:AbstractPrescription,
    C<:AbstractComponentSelection,
} <: AbstractVertex
    quantity::Q
    value::P
    components::C
end


quantity(load::BoundaryLoad, ::Tuple) = load.quantity

# A prescribed value depends on no field, unless it depends on the current position under
# finite strains (under small strains the current and reference positions coincide).
vertex_fields(::BoundaryLoad, ::SmallStrain) = Symbol[]

function vertex_fields(load::BoundaryLoad, assumption::FiniteStrain)
    if depends_on_current_position(load.value)
        return [motion_field(assumption)]
    else
        return Symbol[]
    end
end

is_affine(load::BoundaryLoad, assumption::AbstractKinematicAssumption) =
    isempty(vertex_fields(load, assumption))

subvertices(::BoundaryLoad) = ()


"""
    RobinFlux(field, coefficient, reference, components, configuration)

The flux ``h (T - T_∞)`` of a Robin condition on `field`. Its tag is a normal flux for a
scalar field and a traction for a vector field, per unit area of `configuration`.

# Fields
- `field::Symbol`: the field.
- `coefficient::Float64`: the coefficient ``h``.
- `reference::P`: the reference value ``T_∞``, an [`AbstractPrescription`](@ref).
- `components::C`: the [`AbstractComponentSelection`](@ref).
- `configuration::K`: the configuration whose area the coefficient is measured per unit
  of.
"""
struct RobinFlux{
    P<:AbstractPrescription,
    C<:AbstractComponentSelection,
    K<:AbstractConfiguration,
} <: AbstractVertex
    field::Symbol
    coefficient::Float64
    reference::P
    components::C
    configuration::K
end


function quantity(flux::RobinFlux, fields::Tuple)
    # The tag follows the field's shape and the configuration of the coefficient.
    field = lookup_field(fields, flux.field)
    if field === nothing
        return nothing
    end
    return robin_quantity(field_shape(field), flux.configuration)
end


"""
    robin_quantity(field_shape, configuration)

Return the tag of the flux of a Robin condition on a field of shape `field_shape`, per
unit area of `configuration`, or `nothing` for shapes without one.
"""
robin_quantity(::AbstractFieldShape, ::AbstractConfiguration) = nothing
robin_quantity(::ScalarField, ::ReferenceConfiguration) = ReferenceNormalFlux()
robin_quantity(::ScalarField, ::CurrentConfiguration) = SpatialNormalFlux()
robin_quantity(::VectorField, ::ReferenceConfiguration) = ReferenceTraction()
robin_quantity(::VectorField, ::CurrentConfiguration) = SpatialTraction()

# The flux is affine in its field, with a constant coefficient.
vertex_fields(flux::RobinFlux, ::AbstractKinematicAssumption) = [flux.field]
is_affine(::RobinFlux, ::AbstractKinematicAssumption) = true
subvertices(flux::RobinFlux) = (Value(flux.field),)
