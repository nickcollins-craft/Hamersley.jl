#=
Compressible neo-Hookean hyperelasticity under finite strains.

This is the declaration of the first finite-strain law in the library: its interface (a
deformation gradient in, a first Piola–Kirchhoff stress out) and its symmetric tangent. The
evaluation itself is not yet written.
=#


@doc raw"""
    NeoHookean(; shear_modulus, first_lame_parameter)

A compressible neo-Hookean material, with stored energy

```math
Ψ(F) = \frac{μ}{2} (\operatorname{tr}(F^T F) - 3) - μ \ln J + \frac{λ}{2} (\ln J)^2,
\qquad J = \det F,
```

and first Piola–Kirchhoff stress ``P = ∂Ψ/∂F``. The law is nonlinear in ``F`` and its
tangent ``∂^2 Ψ / ∂F ∂F`` is symmetric.

# Fields
- `μ::T`: the shear modulus (Pa).
- `λ::T`: the first Lamé parameter (Pa).

# References
Holzapfel, G. A. (2000), *Nonlinear Solid Mechanics*, Chapter 6.
"""
struct NeoHookean{T<:Real} <: AbstractLaw
    μ::T
    λ::T
end


function NeoHookean(; shear_modulus::Real, first_lame_parameter::Real)
    # We promote the two parameters to a common type so that the struct is concrete.
    μ, λ = promote(shear_modulus, first_lame_parameter)
    law = NeoHookean(μ, λ)
    return law
end


law_inputs(::NeoHookean) = (F = DeformationGradientTensor(),)
law_outputs(::NeoHookean) = (P = FirstPiolaStress(),)
tangent_is_symmetric(::NeoHookean) = true
