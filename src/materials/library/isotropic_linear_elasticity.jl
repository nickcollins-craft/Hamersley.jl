#=
Isotropic linear elasticity under small strains.

This is the declaration of the simplest law in the library: its interface (an
infinitesimal strain in, a Cauchy stress out), its linearity and its potential. The
evaluation itself is not yet written.
=#


@doc raw"""
    IsotropicLinearElasticity(; first_lame_parameter, shear_modulus)

Isotropic linear elasticity,

```math
σ = λ \operatorname{tr}(ε) I + 2 μ ε,
```

which derives from the stored energy
``Ψ = \tfrac{λ}{2} (\operatorname{tr} ε)^2 + μ \, ε : ε``.
The law is affine in ``ε``, so a problem built only from it (under small strains) is
linear.

# Fields
- `λ::T`: the first Lamé parameter (Pa).
- `μ::T`: the shear modulus (Pa).
"""
struct IsotropicLinearElasticity{T<:Real} <: AbstractLaw
    λ::T
    μ::T
end


function IsotropicLinearElasticity(; first_lame_parameter::Real, shear_modulus::Real)
    # We promote the two parameters to a common type so that the struct is concrete.
    λ, μ = promote(first_lame_parameter, shear_modulus)
    law = IsotropicLinearElasticity(λ, μ)
    return law
end


law_inputs(::IsotropicLinearElasticity) = (ε = InfinitesimalStrainTensor(),)
law_outputs(::IsotropicLinearElasticity) = (σ = CauchyStress(),)
is_affine(::IsotropicLinearElasticity) = true
has_potential(::IsotropicLinearElasticity) = true
