#=
Kinematic assumptions.

A Hamersley problem declares, once and explicitly, whether it is posed under small
strains or finite strains. Everything that depends on the distinction between the
reference and current configurations (conjugate pairs, pull-backs, which measure a stress
can be integrated over) consults this assumption rather than guessing.

Under small strains the deformation gradient is ``F ≈ I`` and ``J ≈ 1``, so the reference
and current configurations coincide to first order, and all stress measures coincide with
the Cauchy stress. Under finite strains they are distinct, and we solve the problem in a
total Lagrangian manner: every integral is taken over the reference configuration
``Ω_0``. "Total Lagrangian" is therefore *how* finite strains are solved, not the name of
the assumption.
=#


"""
    AbstractKinematicAssumption

Supertype of the kinematic assumptions a problem can be posed under.
"""
abstract type AbstractKinematicAssumption end


"""
    SmallStrain()

The small-strain (geometrically linear) assumption. The reference and current
configurations are identified, every pull-back and push-forward is the identity, and the
infinitesimal strain ``ε = \\operatorname{sym}(∇u)`` is work-conjugate to the Cauchy stress
``σ``.
"""
struct SmallStrain <: AbstractKinematicAssumption end


"""
    FiniteStrain(; motion)

The finite-strain (geometrically nonlinear) assumption, solved in a total Lagrangian
manner.

# Fields
- `motion::Symbol`: the name of the displacement field ``u`` that defines the motion
  ``x = X + u(X)``, and hence ``F = I + ∇_0 u`` and ``J = \\det F``. Every pull-back and
  push-forward in the problem uses this motion unless it names another one.
"""
struct FiniteStrain <: AbstractKinematicAssumption
    motion::Symbol
end

FiniteStrain(; motion::Symbol) = FiniteStrain(motion)


"""
    motion_field(assumption)

Return the name of the field that defines the motion, or `nothing` under small strains
(where the motion does not enter any pull-back).
"""
motion_field(::SmallStrain) = nothing
motion_field(assumption::FiniteStrain) = assumption.motion


"""
    effective_configuration(assumption, configuration_value)

Return the configuration that `configuration_value` is equivalent to under `assumption`.

Under [`SmallStrain`](@ref) every configuration except
[`ConfigurationIndependent`](@ref) collapses onto the reference configuration (because
``F ≈ I``). Under [`FiniteStrain`](@ref) configurations are kept distinct.
"""
function effective_configuration(::SmallStrain, configuration_value::AbstractConfiguration)
    # Configuration-independent quantities stay so; everything else is identified with Ω_0.
    if configuration_value isa ConfigurationIndependent
        return configuration_value
    else
        return ReferenceConfiguration()
    end
end

function effective_configuration(
    ::FiniteStrain,
    configuration_value::AbstractConfiguration,
)
    return configuration_value
end


"""
    same_configuration(assumption, a, b)

Return `true` when `a` and `b` (configurations, measures or quantity tags) belong to the
same configuration once `assumption` has been taken into account.
"""
function same_configuration(assumption::AbstractKinematicAssumption, a, b)
    # We first reduce both arguments to configurations, then compare their effective forms.
    configuration_a = effective_configuration(assumption, configuration(a))
    configuration_b = effective_configuration(assumption, configuration(b))
    return same_configuration(configuration_a, configuration_b)
end
