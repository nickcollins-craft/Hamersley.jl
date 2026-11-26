#=
Boundary traces.

A trace is the restriction of a field (or of one of its derivatives) to the boundary.
Which traces appear in a problem is not declared by hand: integrating the weak form by
parts leaves one boundary term per trace, and which traces those are depends only on the
operators that act on the test function (see `formulation/trace_derivation.jl`). This
file defines the traces themselves and pairs each boundary quantity with the trace it is
conjugate to.
=#


"""
    AbstractTrace

Supertype of the boundary traces of a field.
"""
abstract type AbstractTrace end


"""
    ValueTrace()

The value of the field on the boundary, ``u|_{∂Ω}``. This is the ordinary trace of an
``H^1`` field, conjugate to tractions and normal fluxes.
"""
struct ValueTrace <: AbstractTrace end


"""
    NormalGradientTrace()

The normal derivative of the field on the boundary, ``∂_{n_0} u = (∇_0 u) n_0``, which
appears in second-gradient continua and is conjugate to the double force.
"""
struct NormalGradientTrace <: AbstractTrace end


"""
    NormalComponentTrace()

The normal component of a vector field on the boundary, ``q \\cdot n_0``, which is the
natural trace of an ``H(\\mathrm{div})`` field.
"""
struct NormalComponentTrace <: AbstractTrace end


"""
    TangentialComponentTrace()

The tangential component of a vector field on the boundary, ``q × n_0``, which is the
natural trace of an ``H(\\mathrm{curl})`` field.
"""
struct TangentialComponentTrace <: AbstractTrace end


"""
    EdgeValueTrace()

The value of the field on a sharp edge of the boundary (a corner in two dimensions),
conjugate to the edge forces of second-gradient continua.
"""
struct EdgeValueTrace <: AbstractTrace end


"""
    conjugate_trace(boundary_quantity)

Return the trace that the boundary quantity `boundary_quantity` does work against. This
is how a natural boundary condition infers its trace, so that users only name the
quantity they prescribe.
"""
conjugate_trace(::ReferenceTraction) = ValueTrace()
conjugate_trace(::SpatialTraction) = ValueTrace()
conjugate_trace(::ReferenceNormalFlux) = ValueTrace()
conjugate_trace(::SpatialNormalFlux) = ValueTrace()
conjugate_trace(::DoubleForce) = NormalGradientTrace()
conjugate_trace(::EdgeForce) = EdgeValueTrace()
conjugate_trace(::BoundaryPressure) = NormalComponentTrace()


"""
    is_standard_trace(trace)

Return `true` for the ordinary value trace, and `false` for the non-standard traces of
higher-order and mixed formulations. Strict boundary checking in the formulation verification insists that
non-standard traces are dealt with explicitly.
"""
is_standard_trace(::AbstractTrace) = false
is_standard_trace(::ValueTrace) = true
