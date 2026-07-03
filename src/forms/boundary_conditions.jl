abstract type AbstractBoundaryCondition end

struct EssentialBC{T,V} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    boundary::Symbol
    value::V
end

struct NaturalBC{T,Q} <: AbstractBoundaryCondition
    field::Symbol
    trace::T
    boundary::Symbol
    quantity::Q
end

function EssentialBC(
    field::Symbol,
    trace::AbstractTrace;
    on::Symbol,
    value,
)
    return EssentialBC{typeof(trace),typeof(value)}(field, trace, on, value)
end

function NaturalBC(
    field::Symbol,
    trace::AbstractTrace;
    on::Symbol,
    quantity::AbstractBoundaryQuantity,
)
    return NaturalBC{typeof(trace),typeof(quantity)}(field, trace, on, quantity)
end

function NaturalBC(
    field::Symbol;
    on::Symbol,
    quantity::AbstractBoundaryQuantity,
)
    return NaturalBC{Nothing,typeof(quantity)}(field, nothing, on, quantity)
end

bc_field(bc::AbstractBoundaryCondition) = bc.field
bc_trace(bc::AbstractBoundaryCondition) = bc.trace
bc_boundary(bc::AbstractBoundaryCondition) = bc.boundary
bc_quantity(bc::NaturalBC) = bc.quantity
bc_value(bc::EssentialBC) = bc.value
