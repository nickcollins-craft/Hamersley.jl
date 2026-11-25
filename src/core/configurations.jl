#=
Configurations and integration measures.

Every quantity in a Hamersley formulation is "per unit something" of a particular
configuration of the body: the Cauchy stress is a force per unit *current* area, the first
Piola–Kirchhoff stress a force per unit *reference* area, and so on. This file defines
the configurations themselves, and the measures (volume, boundary, interface) that weak
terms are integrated against. The kinematic assumption of a problem (see
`kinematic_assumptions.jl`) decides which of these distinctions actually matter.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Configurations ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractConfiguration

Supertype of the configurations of the body in which a quantity can be expressed.
"""
abstract type AbstractConfiguration end


"""
    ReferenceConfiguration()

The reference (undeformed) configuration ``Ω_0``, with material coordinates ``X``.
"""
struct ReferenceConfiguration <: AbstractConfiguration end


"""
    CurrentConfiguration()

The current (deformed) configuration ``Ω``, with spatial coordinates ``x = φ(X, t)``.
"""
struct CurrentConfiguration <: AbstractConfiguration end


"""
    IntermediateConfiguration()

An intermediate (stress-free) configuration, such as the one produced by the plastic part
``F_p`` of a multiplicative split ``F = F_e F_p``.
"""
struct IntermediateConfiguration <: AbstractConfiguration end


"""
    ConfigurationIndependent()

Marker for quantities that do not belong to any configuration, such as the temperature
``θ`` or the volume ratio ``J``. A configuration-independent quantity is compatible with
every configuration.
"""
struct ConfigurationIndependent <: AbstractConfiguration end


# ---------------------------------------------------------------------------------------- #
# --------------------------------------- Measures --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractMeasure

Supertype of the measures that weak terms are integrated against.
"""
abstract type AbstractMeasure end


"""
    AbstractVolumeMeasure

Measures over the interior of the body (cells of the mesh).
"""
abstract type AbstractVolumeMeasure <: AbstractMeasure end


"""
    AbstractBoundaryMeasure

Measures over the external boundary of the body (boundary facets of the mesh).
"""
abstract type AbstractBoundaryMeasure <: AbstractMeasure end


"""
    AbstractInterfaceMeasure

Measures over interior facets shared by two cells, as needed by discontinuous Galerkin and
interior-penalty methods.
"""
abstract type AbstractInterfaceMeasure <: AbstractMeasure end


"""
    ReferenceVolume()

The reference volume element ``\\mathrm{d}V_0``.
"""
struct ReferenceVolume <: AbstractVolumeMeasure end


"""
    CurrentVolume()

The current volume element ``\\mathrm{d}v = J \\, \\mathrm{d}V_0``.
"""
struct CurrentVolume <: AbstractVolumeMeasure end


"""
    ReferenceBoundary()

The reference area element ``\\mathrm{d}A_0``.
"""
struct ReferenceBoundary <: AbstractBoundaryMeasure end


"""
    CurrentBoundary()

The current area element ``\\mathrm{d}a``, related to the reference one by Nanson's formula
``n \\, \\mathrm{d}a = J F^{-T} N \\, \\mathrm{d}A_0``.
"""
struct CurrentBoundary <: AbstractBoundaryMeasure end


"""
    ReferenceInterface()

The reference area element on interior facets. This is a placeholder: interface terms are
part of the structure so that discontinuous methods can be added later, but they cannot
yet be assembled.
"""
struct ReferenceInterface <: AbstractInterfaceMeasure end


# ---------------------------------------------------------------------------------------- #
# ------------------------------ Configuration bookkeeping ------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    configuration(x)

Return the configuration that `x` belongs to. For a configuration this is `x` itself; for a
measure it is the configuration whose volume or area the measure counts.
"""
configuration(configuration_value::AbstractConfiguration) = configuration_value
configuration(::ReferenceVolume) = ReferenceConfiguration()
configuration(::CurrentVolume) = CurrentConfiguration()
configuration(::ReferenceBoundary) = ReferenceConfiguration()
configuration(::CurrentBoundary) = CurrentConfiguration()
configuration(::ReferenceInterface) = ReferenceConfiguration()


"""
    same_configuration(a::AbstractConfiguration, b::AbstractConfiguration)

Return `true` when `a` and `b` are the same configuration, or when either of them is
[`ConfigurationIndependent`](@ref).

This comparison knows nothing about the kinematic assumption of a problem; the
assumption-aware version, which is what the audit uses, is
`same_configuration(assumption, a, b)` in `kinematic_assumptions.jl`.
"""
function same_configuration(a::AbstractConfiguration, b::AbstractConfiguration)
    # A configuration-independent quantity is compatible with anything.
    if a isa ConfigurationIndependent || b isa ConfigurationIndependent
        return true
    else
        return typeof(a) === typeof(b)
    end
end
