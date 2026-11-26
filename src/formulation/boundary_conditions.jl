#=
Boundary conditions.

An essential condition prescribes a trace of a field; a natural condition prescribes the
boundary quantity that does work against a trace. Which traces exist for a field is not
declared: it is derived from the weak form (see `trace_derivation.jl`). A natural
condition infers its trace from the quantity it prescribes, so users only write what they
are prescribing.

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


"""
    EssentialBC(field; on, value, trace = ValueTrace())

Prescribe the trace `trace` of `field` on the boundary set `on`.

# Fields
- `field::Symbol`: the field constrained.
- `trace::T`: the [`AbstractTrace`](@ref) prescribed.
- `boundary::Symbol`: the name of the facet set.
- `value::V`: the prescribed value (for now often a [`PlaceholderValue`](@ref)).
"""
struct EssentialBC{T<:AbstractTrace,V} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    boundary::Symbol
    value::V
end


function EssentialBC(field::Symbol; on::Symbol, value, trace::AbstractTrace = ValueTrace())
    bc = EssentialBC(field, trace, on, value)
    return bc
end


"""
    NaturalBC(field; on, quantity, trace = conjugate_trace(quantity))

Prescribe the boundary quantity `quantity` on the boundary set `on`, doing work against
a trace of `field`. The trace is inferred from the quantity (a traction works against the
value, a double force against the normal derivative), and may be given explicitly only to
document intent.

# Fields
- `field::Symbol`: the field whose test function the load does work against.
- `trace::T`: the [`AbstractTrace`](@ref) the quantity is conjugate to.
- `boundary::Symbol`: the name of the facet set.
- `quantity::Q`: the [`AbstractBoundaryQuantity`](@ref) prescribed.
"""
struct NaturalBC{T<:AbstractTrace,Q<:AbstractBoundaryQuantity} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    boundary::Symbol
    quantity::Q
end


function NaturalBC(
    field::Symbol;
    on::Symbol,
    quantity::AbstractBoundaryQuantity,
    trace::AbstractTrace = conjugate_trace(quantity),
)
    bc = NaturalBC(field, trace, on, quantity)
    return bc
end


bc_field(bc::AbstractBoundaryCondition) = bc.field
bc_trace(bc::AbstractBoundaryCondition) = bc.trace
bc_boundary(bc::AbstractBoundaryCondition) = bc.boundary
bc_quantity(bc::NaturalBC) = bc.quantity
bc_value(bc::EssentialBC) = bc.value


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

The condition is imposed weakly, by Nitsche's method or a penalty. Not yet implemented.
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
imposition(::NaturalBC, ::Field) = NaturalImposition()

function imposition(bc::EssentialBC, field::Field)
    # The interpolation decides whether the trace is a degree of freedom.
    conformity = family_conformity(field_interpolation(field))
    if can_impose_strongly(conformity, bc.trace)
        return StrongImposition()
    else
        return WeakImposition()
    end
end
