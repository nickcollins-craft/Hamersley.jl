#=
Weak terms.

A weak term pairs a test-side quantity with a flux, integrated against a measure:

```math
\int δ(\text{test}) : \text{flux} \, \mathrm{d}(\text{measure}).
```

The test side is always written explicitly as a `Variation`, so that, for example,
``∫ P : ∇_0 δu \, \mathrm{d}V_0`` and ``∫ S : δE \, \mathrm{d}V_0`` are different terms
with different linearisations. The kind of term (volume, boundary or interface) follows
from its measure.
=#


"""
    AbstractWeakTerm

Supertype of the terms of a weak form.
"""
abstract type AbstractWeakTerm end


"""
    Term(test, flux, measure; on = nothing)

The weak term ``∫ δ(\\text{test}) : \\text{flux} \\, \\mathrm{d}(\\text{measure})``.

# Fields
- `test::T`: the test side, a [`Variation`](@ref).
- `flux::F`: the node paired with the test side.
- `measure::M`: the measure integrated against; it decides whether this is a volume,
  boundary or interface term.
- `domain::Union{Symbol,Nothing}`: the cell set (for volume terms) or facet set (for
  boundary terms) the term is restricted to, or `nothing` for the whole domain or
  boundary.

# Examples
```julia
F = DeformationGradient(:u)
Term(Variation(F), output(mechanics, :P), ReferenceVolume())
```
"""
struct Term{T<:Variation,F<:AbstractNode,M<:AbstractMeasure} <: AbstractWeakTerm
    test::T
    flux::F
    measure::M
    domain::Union{Symbol,Nothing}
end


Term(
    test::Variation,
    flux::AbstractNode,
    measure::AbstractMeasure;
    on::Union{Symbol,Nothing} = nothing,
) = Term(test, flux, measure, on)


term_test(term::Term) = term.test
term_flux(term::Term) = term.flux
term_measure(term::Term) = term.measure
term_domain(term::Term) = term.domain


is_volume_term(term::Term) = term.measure isa AbstractVolumeMeasure
is_boundary_term(term::Term) = term.measure isa AbstractBoundaryMeasure
is_interface_term(term::Term) = term.measure isa AbstractInterfaceMeasure
