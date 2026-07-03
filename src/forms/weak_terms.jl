abstract type AbstractWeakTerm end

struct TestValueTerm{R<:AbstractResponse} <: AbstractWeakTerm
    response::R
end

struct TestGradientTerm{R<:AbstractResponse} <: AbstractWeakTerm
    response::R
end

struct BoundaryTraceTerm{T<:AbstractTrace,Q<:AbstractBoundaryQuantity} <: AbstractWeakTerm
    trace::T
    quantity::Q
end

struct TimeDerivativeTerm{R<:StorageResponse} <: AbstractWeakTerm
    response::R
end

term_response(term::Union{TestValueTerm,TestGradientTerm,TimeDerivativeTerm}) = term.response
term_trace(term::BoundaryTraceTerm) = term.trace
term_quantity(term::BoundaryTraceTerm) = term.quantity
