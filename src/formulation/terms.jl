#=
Weak terms.

A weak term pairs a test-side quantity with a flux, integrated against a measure:

```math
\int δ(\text{test}) : \text{flux} \, \mathrm{d}(\text{measure}).
```

The test side is always written explicitly as a `Variation`, so that, for example,
``∫ P : ∇_0 δu \, \mathrm{d}v_0`` and ``∫ S : δE \, \mathrm{d}v_0`` are different terms
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
- `flux::F`: the vertex paired with the test side.
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
struct Term{T<:Variation,F<:AbstractVertex,M<:AbstractMeasure} <: AbstractWeakTerm
    test::T
    flux::F
    measure::M
    domain::Union{Symbol,Nothing}
end


Term(
    test::Variation,
    flux::AbstractVertex,
    measure::AbstractMeasure;
    on::Union{Symbol,Nothing}=nothing,
) = Term(test, flux, measure, on)


term_test(term::Term) = term.test
term_flux(term::Term) = term.flux
term_measure(term::Term) = term.measure
term_domain(term::Term) = term.domain


is_volume_term(term::Term) = term.measure isa AbstractVolumeMeasure
is_boundary_term(term::Term) = term.measure isa AbstractBoundaryMeasure
is_interface_term(term::Term) = term.measure isa AbstractInterfaceMeasure


"""
    is_affine(term, assumption)

Return `true` when `term` contributes an affine function of the unknowns to the residual
under `assumption`: its test side must be affine (so that its variation is constant), its
flux must be affine, and its measure must not depend on the unknowns.
"""
function is_affine(term::Term, assumption::AbstractKinematicAssumption)
    # The measure is part of the term: dv = J dv_0 under finite strains.
    test_is_affine = is_affine(term_test(term), assumption)
    flux_is_affine = is_affine(term_flux(term), assumption)
    measure_is_affine = is_affine(term_measure(term), assumption)
    return test_is_affine && flux_is_affine && measure_is_affine
end


"""
    term_dependencies(term, assumption)

Return a vector of the names of the fields that the contribution of `term` to the residual
depends on under `assumption`: those of its flux, those of its measure, and those of its
test side when its variation depends on the solution
(``δE = \\operatorname{sym}(F^T ∇_0 δu)`` does; ``δp`` and ``δF = ∇_0 δu`` do not).
"""
function term_dependencies(term::Term, assumption::AbstractKinematicAssumption)
    # The flux and the measure always enter the integrand.
    dependencies = copy(vertex_fields(term_flux(term), assumption))
    for field in measure_fields(term_measure(term), assumption)
        push_unique!(dependencies, field)
    end

    # The test side only adds dependencies when its variation is not constant.
    test = term_test(term)
    if !is_affine(test, assumption)
        for field in vertex_fields(test, assumption)
            push_unique!(dependencies, field)
        end
    end
    return dependencies
end
