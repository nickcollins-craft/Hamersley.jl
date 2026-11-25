#=
Quantity tags.

A quantity tag records what a quantity *means*: which physical object it is, which
configuration it is measured per unit of, and how it transforms between configurations.
It says nothing about where the value comes from. The same tag (say `CauchyStress()`) is
carried by the output of a constitutive law, by an interpolated stress field in a mixed
formulation, or by a pulled-back quantity, and the audit checks conjugacy and
configurations on tags alone.

This file defines the tags, their configurations, the pull-back and push-forward tables,
how tags collapse under small strains, and which pairs of tags are work-conjugate.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Tag hierarchy ------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractQuantity

Supertype of every quantity tag.
"""
abstract type AbstractQuantity end


"""
    AbstractFieldValue{S}

Tags for the point value of a field of shape `S`. Fields that declare a physical meaning
(`represents = Temperature()`) carry that tag; fields that do not carry
[`FieldValue{S}`](@ref FieldValue).
"""
abstract type AbstractFieldValue{S<:AbstractFieldShape} <: AbstractQuantity end


"""
    AbstractKinematicQuantity

Tags for kinematic quantities derived from fields by differential operators.
"""
abstract type AbstractKinematicQuantity <: AbstractQuantity end


"""
    AbstractStressMeasure

Tags for stress measures (and the higher-order stresses of gradient theories).
"""
abstract type AbstractStressMeasure <: AbstractQuantity end


"""
    AbstractFluxQuantity

Tags for vector fluxes of scalar quantities, such as heat or pore-fluid fluxes.
"""
abstract type AbstractFluxQuantity <: AbstractQuantity end


"""
    AbstractDensityQuantity

Tags for quantities per unit volume, such as body forces, heat sources and stored
contents.
"""
abstract type AbstractDensityQuantity <: AbstractQuantity end


"""
    AbstractBoundaryQuantity

Tags for quantities prescribed on the boundary by natural boundary conditions.
"""
abstract type AbstractBoundaryQuantity <: AbstractQuantity end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Field values -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    FieldValue{S}()

The point value of a field of shape `S` that does not declare a more specific meaning.
"""
struct FieldValue{S} <: AbstractFieldValue{S} end

FieldValue(::S) where {S<:AbstractFieldShape} = FieldValue{S}()


"""
    Temperature()

The absolute temperature ``θ``.
"""
struct Temperature <: AbstractFieldValue{ScalarField} end


"""
    Pressure()

A pressure ``p``: a pore pressure, or the hydrostatic pressure of a mixed
displacement–pressure formulation.
"""
struct Pressure <: AbstractFieldValue{ScalarField} end


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Kinematic quantities -------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    MaterialGradient{S}()

The gradient ``∇_0(\\cdot)`` with respect to the material coordinates ``X`` of a field of
shape `S`.
"""
struct MaterialGradient{S} <: AbstractKinematicQuantity end

MaterialGradient(::S) where {S<:AbstractFieldShape} = MaterialGradient{S}()


"""
    SpatialGradient{S}()

The gradient ``∇_x(\\cdot)`` with respect to the spatial coordinates ``x`` of a field of
shape `S`.
"""
struct SpatialGradient{S} <: AbstractKinematicQuantity end

SpatialGradient(::S) where {S<:AbstractFieldShape} = SpatialGradient{S}()


"""
    MaterialDivergence()

The divergence ``\\operatorname{Div}(\\cdot)`` of a vector field with respect to the
material coordinates.
"""
struct MaterialDivergence <: AbstractKinematicQuantity end


"""
    MaterialSecondGradient()

The second gradient ``∇_0 ∇_0 u`` of a vector field, as it appears in second-gradient
continua.
"""
struct MaterialSecondGradient <: AbstractKinematicQuantity end


"""
    DeformationGradientTensor()

The deformation gradient ``F = I + ∇_0 u``, a two-point tensor.
"""
struct DeformationGradientTensor <: AbstractKinematicQuantity end


"""
    GreenLagrangeStrainTensor()

The Green–Lagrange strain ``E = \\tfrac{1}{2}(F^T F - I)``, a material tensor.
"""
struct GreenLagrangeStrainTensor <: AbstractKinematicQuantity end


"""
    InfinitesimalStrainTensor()

The infinitesimal strain ``ε = \\operatorname{sym}(∇u)``. It is only meaningful under
[`SmallStrain`](@ref).
"""
struct InfinitesimalStrainTensor <: AbstractKinematicQuantity end


"""
    JacobianDeterminant()

The volume ratio ``J = \\det F = \\mathrm{d}v / \\mathrm{d}V_0``.
"""
struct JacobianDeterminant <: AbstractKinematicQuantity end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Stress measures ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    CauchyStress()

The Cauchy (true) stress ``σ``: force per unit current area.
"""
struct CauchyStress <: AbstractStressMeasure end


"""
    KirchhoffStress()

The Kirchhoff stress ``τ = J σ``: a spatial tensor, but per unit *reference* volume.
"""
struct KirchhoffStress <: AbstractStressMeasure end


"""
    FirstPiolaStress()

The first Piola–Kirchhoff stress ``P = J σ F^{-T}``: force per unit reference area, a
two-point tensor.
"""
struct FirstPiolaStress <: AbstractStressMeasure end


"""
    SecondPiolaStress()

The second Piola–Kirchhoff stress ``S = F^{-1} P``, a material tensor.
"""
struct SecondPiolaStress <: AbstractStressMeasure end


"""
    HyperStress()

The hyperstress (double stress) of a second-gradient continuum, conjugate to the second
gradient ``∇_0 ∇_0 u``.
"""
struct HyperStress <: AbstractStressMeasure end


# ---------------------------------------------------------------------------------------- #
# --------------------------------- Fluxes and densities --------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    MaterialFlux()

A flux ``Q`` per unit reference area, such as the material heat flux or the material
Darcy flux.
"""
struct MaterialFlux <: AbstractFluxQuantity end


"""
    SpatialFlux()

A flux ``q`` per unit current area, related to the material flux by the Piola transform
``Q = J F^{-1} q``.
"""
struct SpatialFlux <: AbstractFluxQuantity end


"""
    ReferenceVolumeDensity{S}()

A quantity of shape `S` per unit reference volume, such as a body force ``b_0``, a heat
source ``r_0`` or a stored content.
"""
struct ReferenceVolumeDensity{S} <: AbstractDensityQuantity end

ReferenceVolumeDensity(::S) where {S<:AbstractFieldShape} = ReferenceVolumeDensity{S}()


"""
    SpatialVolumeDensity{S}()

A quantity of shape `S` per unit current volume, related to the reference density by
``r_0 = J r``.
"""
struct SpatialVolumeDensity{S} <: AbstractDensityQuantity end

SpatialVolumeDensity(::S) where {S<:AbstractFieldShape} = SpatialVolumeDensity{S}()


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Boundary quantities --------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    ReferenceTraction()

A traction ``T = P N`` per unit reference area.
"""
struct ReferenceTraction <: AbstractBoundaryQuantity end


"""
    SpatialTraction()

A traction ``t = σ n`` per unit current area. Under finite strains this is a follower
load, pulled back with Nanson's formula.
"""
struct SpatialTraction <: AbstractBoundaryQuantity end


"""
    DoubleForce()

The double force of a second-gradient continuum, conjugate to the normal derivative
``∂_N u`` on the boundary.
"""
struct DoubleForce <: AbstractBoundaryQuantity end


"""
    EdgeForce()

The line (wedge) force that second-gradient continua carry on sharp edges of the
boundary.
"""
struct EdgeForce <: AbstractBoundaryQuantity end


"""
    ReferenceNormalFlux()

A scalar normal flux ``Q \\cdot N`` per unit reference area.
"""
struct ReferenceNormalFlux <: AbstractBoundaryQuantity end


"""
    SpatialNormalFlux()

A scalar normal flux ``q \\cdot n`` per unit current area.
"""
struct SpatialNormalFlux <: AbstractBoundaryQuantity end


"""
    BoundaryPressure()

A prescribed pressure on the boundary. In flux-based mixed formulations (such as mixed
Darcy flow) the pressure is the *natural* condition, conjugate to the normal component of
the flux test function.
"""
struct BoundaryPressure <: AbstractBoundaryQuantity end


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Time rates -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    RateOf{Q}(quantity)

The time rate ``\\dot{q}`` of a quantity tagged `Q`. A rate has the same configuration and
the same conjugates as the quantity itself; how it is computed from the unknowns is the
business of the time integrator.

# Fields
- `quantity::Q`: the tag of the quantity whose rate is taken.
"""
struct RateOf{Q<:AbstractQuantity} <: AbstractQuantity
    quantity::Q
end


"""
    strip_rate(tag)

Return the tag underneath any number of [`RateOf`](@ref) wrappers.
"""
strip_rate(tag::AbstractQuantity) = tag
strip_rate(tag::RateOf) = strip_rate(tag.quantity)


# ---------------------------------------------------------------------------------------- #
# ----------------------------------- Placeholder values --------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    PlaceholderValue(name)

A named stand-in for a prescribed value (for example the boundary temperature
``\\bar{θ}``) while the representation of prescribed values is not yet settled.

# Fields
- `name::Symbol`: the symbol of the prescribed value.
"""
struct PlaceholderValue
    name::Symbol
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Configurations ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    quantity_configuration(tag)

Return the configuration that the quantity tagged `tag` is measured per unit of. This is
the configuration whose measure the quantity can be integrated against without a
pull-back.

Note that the Kirchhoff stress is a spatial tensor but is measured per unit *reference*
volume, which is why it belongs to the reference configuration here.
"""
quantity_configuration(::AbstractFieldValue) = ConfigurationIndependent()
quantity_configuration(::MaterialGradient) = ReferenceConfiguration()
quantity_configuration(::SpatialGradient) = CurrentConfiguration()
quantity_configuration(::MaterialDivergence) = ReferenceConfiguration()
quantity_configuration(::MaterialSecondGradient) = ReferenceConfiguration()
quantity_configuration(::DeformationGradientTensor) = ReferenceConfiguration()
quantity_configuration(::GreenLagrangeStrainTensor) = ReferenceConfiguration()
quantity_configuration(::InfinitesimalStrainTensor) = ReferenceConfiguration()
quantity_configuration(::JacobianDeterminant) = ConfigurationIndependent()
quantity_configuration(::CauchyStress) = CurrentConfiguration()
quantity_configuration(::KirchhoffStress) = ReferenceConfiguration()
quantity_configuration(::FirstPiolaStress) = ReferenceConfiguration()
quantity_configuration(::SecondPiolaStress) = ReferenceConfiguration()
quantity_configuration(::HyperStress) = ReferenceConfiguration()
quantity_configuration(::MaterialFlux) = ReferenceConfiguration()
quantity_configuration(::SpatialFlux) = CurrentConfiguration()
quantity_configuration(::ReferenceVolumeDensity) = ReferenceConfiguration()
quantity_configuration(::SpatialVolumeDensity) = CurrentConfiguration()
quantity_configuration(::ReferenceTraction) = ReferenceConfiguration()
quantity_configuration(::SpatialTraction) = CurrentConfiguration()
quantity_configuration(::DoubleForce) = ReferenceConfiguration()
quantity_configuration(::EdgeForce) = ReferenceConfiguration()
quantity_configuration(::ReferenceNormalFlux) = ReferenceConfiguration()
quantity_configuration(::SpatialNormalFlux) = CurrentConfiguration()
quantity_configuration(::BoundaryPressure) = ConfigurationIndependent()
quantity_configuration(tag::RateOf) = quantity_configuration(tag.quantity)

# Tags take part in the generic configuration bookkeeping of `configurations.jl`.
configuration(tag::AbstractQuantity) = quantity_configuration(tag)


# ---------------------------------------------------------------------------------------- #
# ---------------------------- Transformations between configurations -------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractTransformationRule

Supertype of the rules by which a quantity is mapped between the current and reference
configurations. The rule is a property of the quantity's tensorial character, so users
write `PullBack(q)` and never name the rule themselves.
"""
abstract type AbstractTransformationRule end


"""
    PiolaVectorRule()

The Piola transform of a flux, ``Q = J F^{-1} q``.
"""
struct PiolaVectorRule <: AbstractTransformationRule end


"""
    PiolaStressRule()

The Piola transform of the Cauchy stress, ``P = J σ F^{-T}``.
"""
struct PiolaStressRule <: AbstractTransformationRule end


"""
    KirchhoffStressRule()

The map from the Kirchhoff stress to the first Piola–Kirchhoff stress, ``P = τ F^{-T}``.
"""
struct KirchhoffStressRule <: AbstractTransformationRule end


"""
    CovariantRule()

The covariant map of a gradient, ``∇_0(\\cdot) = F^T ∇_x(\\cdot)``.
"""
struct CovariantRule <: AbstractTransformationRule end


"""
    VolumeDensityRule()

The map of a volume density, ``r_0 = J r``.
"""
struct VolumeDensityRule <: AbstractTransformationRule end


"""
    NansonRule()

The map of an area density by Nanson's formula, ``n \\, \\mathrm{d}a = J F^{-T} N \\,
\\mathrm{d}A_0``.
"""
struct NansonRule <: AbstractTransformationRule end


"""
    transformation_rule(tag)

Return the rule by which the quantity tagged `tag` is pulled back to (or pushed forward
from) the reference configuration, or `nothing` if the quantity has no such map.
"""
transformation_rule(::AbstractQuantity) = nothing
transformation_rule(::SpatialFlux) = PiolaVectorRule()
transformation_rule(::MaterialFlux) = PiolaVectorRule()
transformation_rule(::CauchyStress) = PiolaStressRule()
transformation_rule(::KirchhoffStress) = KirchhoffStressRule()
transformation_rule(::SpatialGradient) = CovariantRule()
transformation_rule(::MaterialGradient) = CovariantRule()
transformation_rule(::SpatialVolumeDensity) = VolumeDensityRule()
transformation_rule(::ReferenceVolumeDensity) = VolumeDensityRule()
transformation_rule(::SpatialTraction) = NansonRule()
transformation_rule(::SpatialNormalFlux) = NansonRule()


"""
    pulled_back(tag)

Return the tag of the reference-configuration counterpart of the quantity tagged `tag`,
or `nothing` if it has none (for example because it already belongs to the reference
configuration).
"""
pulled_back(::AbstractQuantity) = nothing
pulled_back(::CauchyStress) = FirstPiolaStress()
pulled_back(::KirchhoffStress) = FirstPiolaStress()
pulled_back(::SpatialFlux) = MaterialFlux()
pulled_back(::SpatialGradient{S}) where {S} = MaterialGradient{S}()
pulled_back(::SpatialVolumeDensity{S}) where {S} = ReferenceVolumeDensity{S}()
pulled_back(::SpatialTraction) = ReferenceTraction()
pulled_back(::SpatialNormalFlux) = ReferenceNormalFlux()
pulled_back(tag::AbstractFieldValue) = tag
pulled_back(tag::JacobianDeterminant) = tag
pulled_back(tag::RateOf) = rate_of(pulled_back(tag.quantity))


"""
    pushed_forward(tag)

Return the tag of the current-configuration counterpart of the quantity tagged `tag`, or
`nothing` if it has none.
"""
pushed_forward(::AbstractQuantity) = nothing
pushed_forward(::FirstPiolaStress) = CauchyStress()
pushed_forward(::MaterialFlux) = SpatialFlux()
pushed_forward(::MaterialGradient{S}) where {S} = SpatialGradient{S}()
pushed_forward(::ReferenceVolumeDensity{S}) where {S} = SpatialVolumeDensity{S}()
pushed_forward(::ReferenceTraction) = SpatialTraction()
pushed_forward(::ReferenceNormalFlux) = SpatialNormalFlux()
pushed_forward(tag::AbstractFieldValue) = tag
pushed_forward(tag::JacobianDeterminant) = tag
pushed_forward(tag::RateOf) = rate_of(pushed_forward(tag.quantity))


# A rate of "no counterpart" is still no counterpart.
rate_of(::Nothing) = nothing
rate_of(tag::AbstractQuantity) = RateOf(tag)


# ---------------------------------------------------------------------------------------- #
# ----------------------------- Canonical tags per assumption ---------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    canonical_quantity(assumption, tag)

Return the tag that `tag` is equivalent to under `assumption`.

Under [`FiniteStrain`](@ref) every tag is its own canonical form. Under
[`SmallStrain`](@ref) the distinction between configurations vanishes (``F ≈ I``,
``J ≈ 1``), so every stress measure collapses onto the Cauchy stress, every strain
measure onto the infinitesimal strain, and every spatial quantity onto its material
counterpart.
"""
canonical_quantity(::FiniteStrain, tag::AbstractQuantity) = tag
canonical_quantity(::SmallStrain, tag::AbstractQuantity) = tag
canonical_quantity(::SmallStrain, ::KirchhoffStress) = CauchyStress()
canonical_quantity(::SmallStrain, ::FirstPiolaStress) = CauchyStress()
canonical_quantity(::SmallStrain, ::SecondPiolaStress) = CauchyStress()
canonical_quantity(::SmallStrain, ::GreenLagrangeStrainTensor) =
    InfinitesimalStrainTensor()
canonical_quantity(::SmallStrain, ::DeformationGradientTensor) =
    MaterialGradient{VectorField}()
canonical_quantity(::SmallStrain, ::SpatialGradient{S}) where {S} = MaterialGradient{S}()
canonical_quantity(::SmallStrain, ::SpatialFlux) = MaterialFlux()
canonical_quantity(::SmallStrain, ::SpatialVolumeDensity{S}) where {S} =
    ReferenceVolumeDensity{S}()
canonical_quantity(::SmallStrain, ::SpatialTraction) = ReferenceTraction()
canonical_quantity(::SmallStrain, ::SpatialNormalFlux) = ReferenceNormalFlux()
canonical_quantity(assumption::SmallStrain, tag::RateOf) =
    RateOf(canonical_quantity(assumption, tag.quantity))


"""
    requires_small_strain(tag)

Return `true` for quantities that are only meaningful under [`SmallStrain`](@ref), such as
the infinitesimal strain.
"""
requires_small_strain(::AbstractQuantity) = false
requires_small_strain(::InfinitesimalStrainTensor) = true
requires_small_strain(tag::RateOf) = requires_small_strain(tag.quantity)


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Work conjugacy ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    is_conjugate_pair(assumption, a, b)

Return `true` when the canonical tags `a` and `b`, in that order, form a work-conjugate
pair under `assumption`. This is the table that [`are_conjugate`](@ref) consults; it is
written one way round and symmetrised there.

Pairs that hold under either assumption are listed in the two-argument table
`is_conjugate_pair(a, b)`; pairs that hold only under one assumption add a method with
that assumption as the first argument. The table is deliberately conservative: a pair
that is not listed is reported by the audit, and adding a new physical pairing means
adding a method here.
"""
is_conjugate_pair(::AbstractKinematicAssumption, a, b) = is_conjugate_pair(a, b)

# Under small strains ∇δu : σ = ε(δu) : σ, because σ is symmetric.
is_conjugate_pair(::SmallStrain, ::MaterialGradient{VectorField}, ::CauchyStress) = true


# Anything not listed below is not a conjugate pair.
is_conjugate_pair(a, b) = false

# Mechanical pairs: δF : P, ∇_0 δu : P, δE : S, ∇_x δu : σ, ∇_x δu : τ, δε : σ.
is_conjugate_pair(::DeformationGradientTensor, ::FirstPiolaStress) = true
is_conjugate_pair(::MaterialGradient{VectorField}, ::FirstPiolaStress) = true
is_conjugate_pair(::GreenLagrangeStrainTensor, ::SecondPiolaStress) = true
is_conjugate_pair(::SpatialGradient{VectorField}, ::CauchyStress) = true
is_conjugate_pair(::SpatialGradient{VectorField}, ::KirchhoffStress) = true
is_conjugate_pair(::InfinitesimalStrainTensor, ::CauchyStress) = true

# Second-gradient continua: ∇∇δu ⋮ hyperstress.
is_conjugate_pair(::MaterialSecondGradient, ::HyperStress) = true

# Scalar transport: ∇δθ · flux, in either configuration.
is_conjugate_pair(::MaterialGradient{ScalarField}, ::MaterialFlux) = true
is_conjugate_pair(::SpatialGradient{ScalarField}, ::SpatialFlux) = true

# Values against volume densities: δu · b, δθ r, δθ (stored content).
is_conjugate_pair(::AbstractFieldValue{S}, ::ReferenceVolumeDensity{S}) where {S} = true
is_conjugate_pair(::AbstractFieldValue{S}, ::SpatialVolumeDensity{S}) where {S} = true

# Mixed formulations: δp Div u, δp (J - 1), and δq · K⁻¹q (a gradient-like output of a
# resistivity law).
is_conjugate_pair(::AbstractFieldValue{ScalarField}, ::MaterialDivergence) = true
is_conjugate_pair(::AbstractFieldValue{ScalarField}, ::JacobianDeterminant) = true
is_conjugate_pair(::MaterialFlux, ::MaterialGradient{ScalarField}) = true

# Interpolated scalar values against each other, as in penalty or compliance terms δp p/κ.
is_conjugate_pair(::AbstractFieldValue{ScalarField}, ::AbstractFieldValue{ScalarField}) =
    true


"""
    are_conjugate(assumption, test_tag, flux_tag)

Return `true` when the test-side quantity `test_tag` and the flux `flux_tag` of a weak term
form a work-conjugate pair under `assumption`.

Both tags are first reduced to their canonical forms (so that, under small strains,
``ε : σ`` and ``∇u : P`` are recognised as the same pairing), and any time rates are
stripped (a rate has the conjugates of the quantity itself). The table is symmetric, which
is what lets stress-based mixed formulations pair an interpolated stress test function
with a strain.
"""
function are_conjugate(
    assumption::AbstractKinematicAssumption,
    test_tag::AbstractQuantity,
    flux_tag::AbstractQuantity,
)
    # Reduce both tags to the forms the table is written in.
    a = strip_rate(canonical_quantity(assumption, test_tag))
    b = strip_rate(canonical_quantity(assumption, flux_tag))

    # The table is written one way round, so we look the pair up in both orders.
    conjugate = is_conjugate_pair(assumption, a, b) || is_conjugate_pair(assumption, b, a)
    return conjugate
end
