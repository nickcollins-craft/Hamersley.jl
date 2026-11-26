#=
Formulation templates.

Some formulations are discretised differently from how users think about them. A
second-gradient continuum is posed in terms of ``u``, ``∂_{n_0} u`` and the double force,
but with ``C^0`` elements it is discretised as a micromorphic formulation with a relaxed
gradient field ``ψ ≈ ∇u``. The translation between the two belongs to the formulation, not
to each user, so a template declares it once through `translate_trace`.
=#


"""
    AbstractFormulationTemplate

Supertype of formulation templates. A template that discretises a continuum theory
through auxiliary fields overrides [`translate_trace`](@ref).
"""
abstract type AbstractFormulationTemplate end


"""
    NoTemplate()

The default: boundary conditions are written directly on the discretised fields.
"""
struct NoTemplate <: AbstractFormulationTemplate end


"""
    translate_trace(template, field, trace)

Map a boundary condition written in terms of the continuum theory (on `field`, trace
`trace`) onto the discretised formulation. Returns a tuple `(field, trace)`; the default
is the identity.
"""
translate_trace(::AbstractFormulationTemplate, field::Symbol, trace::AbstractTrace) =
    (field, trace)
