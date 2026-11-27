#=
Boundary conditions.

An essential condition prescribes a trace of a field; a natural condition prescribes the
boundary quantity that does work against a trace; a Robin condition loads a trace with a
quantity proportional to it. Which traces exist for a field is not declared: it is derived
from the weak form (see `trace_derivation.jl`). A natural condition infers its trace from
the quantity it prescribes, so users only write what they are prescribing.

Every condition is built from the same parts: *what* is prescribed (the condition's type,
and its trace or quantity), *where* (a region and a selection of components, see
`core/regions.jl`) and *how much* (a value prescription, see `core/prescriptions.jl`).
*When* a condition applies is decided by the problem's schedule (see `schedules.jl`).

Whether an essential condition can be imposed strongly (by constraining degrees of
freedom) depends on whether the field's interpolation can represent the trace, which
follows from Ferrite's conformity of the interpolation. Otherwise it must be imposed
weakly (by Nitsche's method), which is not yet implemented.
=#


"""
    AbstractBoundaryCondition

Supertype of boundary conditions.
"""
abstract type AbstractBoundaryCondition end


# ---------------------------------------------------------------------------------------- #
# --------------------------------- The kinds of condition ------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    EssentialBC(field; on, value, trace = ValueTrace(), components = nothing)

Prescribe the trace `trace` of `field` on the region `on`.

# Fields
- `field::Symbol`: the field constrained.
- `trace::T`: the [`AbstractTrace`](@ref) prescribed.
- `region::R`: the [`AbstractRegion`](@ref); a `Symbol` names a [`FacetSet`](@ref).
- `components::C`: the selected components of the trace (see
  [`component_selection`](@ref)).
- `value::P`: the prescribed value, an [`AbstractPrescription`](@ref); numbers, tuples
  and functions are converted (see [`prescription`](@ref)).
"""
struct EssentialBC{
    T<:AbstractTrace,
    R<:AbstractRegion,
    C<:AbstractComponentSelection,
    P<:AbstractPrescription,
} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    region::R
    components::C
    value::P
end


function EssentialBC(
    field::Symbol;
    on::Union{Symbol,AbstractRegion},
    value,
    trace::AbstractTrace=ValueTrace(),
    components=nothing,
)
    bc = EssentialBC(
        field,
        trace,
        as_region(on),
        component_selection(components),
        prescription(value),
    )
    return bc
end


"""
    NaturalBC(field; on, quantity, value, components = nothing,
              trace = conjugate_trace(quantity))

Prescribe the boundary quantity `quantity` on the region `on`, doing work against a trace
of `field`. The trace is inferred from the quantity (a traction works against the value, a
double force against the normal derivative), and may be given explicitly only to document
intent. The value is measured per unit of the quantity's configuration: a
`ReferenceTraction` per unit reference area, a `SpatialTraction` or `SpatialPressure` per
unit current area (a follower load).

# Fields
- `field::Symbol`: the field whose test function the load does work against.
- `trace::T`: the [`AbstractTrace`](@ref) the quantity is conjugate to.
- `region::R`: the [`AbstractRegion`](@ref); a `Symbol` names a [`FacetSet`](@ref).
- `components::C`: the selected components of the quantity.
- `quantity::Q`: the [`AbstractBoundaryQuantity`](@ref) prescribed.
- `value::P`: the prescribed value, an [`AbstractPrescription`](@ref).
"""
struct NaturalBC{
    T<:AbstractTrace,
    R<:AbstractRegion,
    C<:AbstractComponentSelection,
    Q<:AbstractBoundaryQuantity,
    P<:AbstractPrescription,
} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    region::R
    components::C
    quantity::Q
    value::P
end


function NaturalBC(
    field::Symbol;
    on::Union{Symbol,AbstractRegion},
    quantity::AbstractBoundaryQuantity,
    value,
    components=nothing,
    trace::AbstractTrace=conjugate_trace(quantity),
)
    bc = NaturalBC(
        field,
        trace,
        as_region(on),
        component_selection(components),
        quantity,
        prescription(value),
    )
    return bc
end


"""
    RobinBC(field; on, coefficient, reference = 0.0, components = nothing,
            measure = ReferenceBoundary())

Load the value trace of `field` on the facet region `on` in proportion to its distance
from a reference value: the boundary term ``∫ δT \\, h (T - T_∞) \\, \\mathrm{d}A`` of a
convective face, or ``∫ δu \\cdot k (u - u_∞) \\, \\mathrm{d}A`` of an elastic support.

The coefficient is a constant, because it multiplies the unknown and so enters the
tangent; the reference value is a prescription, because it only enters the right-hand
side. With `measure = CurrentBoundary()` the coefficient is per unit current area, so that
under finite strains the term follows the deformation.

# Fields
- `field::Symbol`: the field loaded.
- `region::R`: the [`AbstractRegion`](@ref); a `Symbol` names a [`FacetSet`](@ref).
- `components::C`: the selected components.
- `coefficient::Float64`: the coefficient ``h``.
- `reference::P`: the reference value ``T_∞``, an [`AbstractPrescription`](@ref).
- `measure::M`: the boundary measure the term is integrated over.
"""
struct RobinBC{
    R<:AbstractRegion,
    C<:AbstractComponentSelection,
    P<:AbstractPrescription,
    M<:AbstractMeasure,
} <: AbstractBoundaryCondition
    field::Symbol
    region::R
    components::C
    coefficient::Float64
    reference::P
    measure::M
end


function RobinBC(
    field::Symbol;
    on::Union{Symbol,AbstractRegion},
    coefficient::Real,
    reference=0.0,
    components=nothing,
    measure::AbstractMeasure=ReferenceBoundary(),
)
    bc = RobinBC(
        field,
        as_region(on),
        component_selection(components),
        Float64(coefficient),
        prescription(reference),
        measure,
    )
    return bc
end


bc_field(bc::AbstractBoundaryCondition) = bc.field
bc_region(bc::AbstractBoundaryCondition) = bc.region
bc_components(bc::AbstractBoundaryCondition) = bc.components
bc_trace(bc::Union{EssentialBC,NaturalBC}) = bc.trace
bc_trace(::RobinBC) = ValueTrace()
bc_quantity(bc::NaturalBC) = bc.quantity
bc_value(bc::Union{EssentialBC,NaturalBC}) = bc.value
bc_value(bc::RobinBC) = bc.reference


"""
    condition_label(bc)

Return a short name for the kind of `bc`, for messages.
"""
condition_label(::EssentialBC) = "an essential condition"
condition_label(::NaturalBC) = "a natural condition"
condition_label(::RobinBC) = "a Robin condition"


"""
    with_components(bc, components)

Return a copy of `bc` acting on the component selection `components` instead. This is how
a condition inherited from an earlier phase is restricted when a later phase replaces
some of its components.
"""
with_components(bc::EssentialBC, components::AbstractComponentSelection) =
    EssentialBC(bc.field, bc.trace, bc.region, components, bc.value)
with_components(bc::NaturalBC, components::AbstractComponentSelection) =
    NaturalBC(bc.field, bc.trace, bc.region, components, bc.quantity, bc.value)
with_components(bc::RobinBC, components::AbstractComponentSelection) =
    RobinBC(bc.field, bc.region, components, bc.coefficient, bc.reference, bc.measure)


# ---------------------------------------------------------------------------------------- #
# --------------------------------- Regions and their traces ----------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    carries_trace(region, trace)

Return `true` when `trace` is one that integration by parts produces on `region`: the
facet traces on facets, and the edge values of second-gradient continua on edges.
Conditions on other regions (a pinned vertex, say) are discrete rather than variational,
and do not enter the status of any trace.
"""
carries_trace(::AbstractRegion, ::AbstractTrace) = false
carries_trace(::FacetSet, ::AbstractTrace) = true
carries_trace(::FacetSet, ::EdgeValueTrace) = false
carries_trace(::EdgeSet, ::EdgeValueTrace) = true


"""
    is_discrete_condition(bc)

Return `true` when `bc` prescribes a value at points or along edges of an ordinary
continuum (an essential value trace on a vertex or edge region), which is not a trace of
the weak form and so is not checked against the derived traces.
"""
is_discrete_condition(bc::AbstractBoundaryCondition) =
    bc isa EssentialBC && bc_trace(bc) isa ValueTrace && !(bc_region(bc) isa FacetSet)


"""
    region_issue(bc)

Return a message if `bc` acts on a kind of region it cannot act on, or `nothing`.
"""
function region_issue(bc::EssentialBC)
    region = bc_region(bc)
    trace = bc_trace(bc)

    # A value can be prescribed anywhere; the edge value of a second-gradient continuum
    # only on edges; every other trace needs a facet, since it involves the normal.
    if trace isa ValueTrace || carries_trace(region, trace)
        return nothing
    end
    return "trace $trace cannot be prescribed on $(region_label(region))"
end

function region_issue(bc::NaturalBC)
    region = bc_region(bc)
    quantity = bc_quantity(bc)

    # Edge forces act on edges; every other boundary quantity is measured per unit area.
    expected = quantity isa EdgeForce ? EdgeSet : FacetSet
    if region isa expected
        return nothing
    elseif region isa VertexSet
        return "$quantity cannot act on $(region_label(region)); concentrated forces " *
               "at points are not yet implemented"
    else
        return "$quantity cannot act on $(region_label(region))"
    end
end

function region_issue(bc::RobinBC)
    # A Robin condition is a boundary integral, so it needs facets.
    region = bc_region(bc)
    if region isa FacetSet
        return nothing
    end
    return "a Robin condition cannot act on $(region_label(region))"
end


# ---------------------------------------------------------------------------------------- #
# --------------------------------------- Value shapes ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    trace_value_shape(field_shape, trace)

Return the shape of the value of `trace` for a field of shape `field_shape`, with every
component selected.
"""
trace_value_shape(::AbstractFieldShape, ::AbstractTrace) = UnknownShape()
trace_value_shape(::ScalarField, ::Union{ValueTrace,NormalGradientTrace,EdgeValueTrace}) =
    ScalarValue()
trace_value_shape(::VectorField, ::Union{ValueTrace,NormalGradientTrace,EdgeValueTrace}) =
    VectorValue(nothing)
trace_value_shape(::AbstractFieldShape, ::NormalComponentTrace) = ScalarValue()
trace_value_shape(::VectorField, ::TangentialComponentTrace) = VectorValue(nothing)


"""
    quantity_value_shape(quantity, field_shape)

Return the shape of the value of the boundary quantity `quantity` acting on a field of
shape `field_shape`, with every component selected.
"""
quantity_value_shape(::AbstractBoundaryQuantity, ::AbstractFieldShape) = UnknownShape()
quantity_value_shape(::Union{ReferenceTraction,SpatialTraction}, ::AbstractFieldShape) =
    VectorValue(nothing)
quantity_value_shape(
    ::Union{
        ReferenceNormalFlux,
        SpatialNormalFlux,
        ReferencePressure,
        SpatialPressure,
        BoundaryPressure,
    },
    ::AbstractFieldShape,
) = ScalarValue()
quantity_value_shape(::Union{DoubleForce,EdgeForce}, field_shape::AbstractFieldShape) =
    trace_value_shape(field_shape, ValueTrace())


"""
    full_value_shape(bc, field_shape)

Return the shape of the value of `bc`, on a field of shape `field_shape`, with every
component selected.
"""
full_value_shape(bc::EssentialBC, field_shape::AbstractFieldShape) =
    trace_value_shape(field_shape, bc_trace(bc))
full_value_shape(bc::NaturalBC, field_shape::AbstractFieldShape) =
    quantity_value_shape(bc_quantity(bc), field_shape)
full_value_shape(::RobinBC, field_shape::AbstractFieldShape) =
    trace_value_shape(field_shape, ValueTrace())


"""
    selected_shape(shape, selection)

Return the shape of a value of shape `shape` restricted to the components `selection`.
"""
selected_shape(shape::AbstractValueShape, ::AbstractComponentSelection) = shape
selected_shape(::VectorValue, selection::AbstractComponentSelection) =
    VectorValue(selected_count(selection))
selected_shape(shape::VectorValue, ::AllComponents) = shape


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Imposition -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractImposition

Supertype of the ways a boundary condition can be imposed.
"""
abstract type AbstractImposition end


"""
    StrongImposition()

The condition is imposed by constraining degrees of freedom (Ferrite's
`ConstraintHandler`).
"""
struct StrongImposition <: AbstractImposition end


"""
    WeakImposition()

The condition is imposed weakly, by Nitsche's method. Not yet implemented.
"""
struct WeakImposition <: AbstractImposition end


"""
    NaturalImposition()

The condition enters the weak form as a boundary integral.
"""
struct NaturalImposition <: AbstractImposition end


"""
    can_impose_strongly(conformity, trace)

Return `true` when an interpolation with Ferrite conformity `conformity` has degrees of
freedom that represent `trace` directly: values for ``H^1``, normal components for
``H(\\mathrm{div})`` and tangential components for ``H(\\mathrm{curl})``.
"""
can_impose_strongly(conformity, trace::AbstractTrace) = false
can_impose_strongly(::Ferrite.H1Conformity, ::ValueTrace) = true
can_impose_strongly(::Ferrite.HdivConformity, ::NormalComponentTrace) = true
can_impose_strongly(::Ferrite.HcurlConformity, ::TangentialComponentTrace) = true


"""
    imposition(bc, field)

Return how `bc` is imposed on `field`.
"""
imposition(::Union{NaturalBC,RobinBC}, ::Field) = NaturalImposition()

function imposition(bc::EssentialBC, field::Field)
    # The interpolation decides whether the trace is a degree of freedom.
    conformity = family_conformity(field_interpolation(field))
    if can_impose_strongly(conformity, bc.trace)
        return StrongImposition()
    else
        return WeakImposition()
    end
end
