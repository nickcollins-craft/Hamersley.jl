#=
Linearity.

A problem is linear when its residual is affine in the unknowns. That property is inferred
from the graph (every term's test side and flux must be affine, and no boundary load may
follow the deformation), so a user's `linear = true` is checked once, in the audit, and
then trusted by the solver.

Linearity is one of two independent properties that let a solver skip work. The other is
that the tangent has not changed since it was last factorised, which is what
`FactorisationKey` records. A linear problem with a fixed time step needs one
factorisation for the whole analysis; a linear problem with adaptive time steps
refactorises whenever the step changes but never iterates; and a nonlinear problem can
reuse a factorisation (modified Newton) without being linear.
=#


"""
    is_affine(problem)

Return `true` when the residual of `problem` is an affine function of its unknowns.
"""
function is_affine(problem::Problem)
    assumption = kinematic_assumption(problem)

    # Every term must have an affine test side (so that its variation is constant) and an
    # affine flux.
    for (_, term) in problem_terms(problem)
        test_is_affine = is_affine(term_test(term), assumption)
        flux_is_affine = is_affine(term_flux(term), assumption)
        if !(test_is_affine && flux_is_affine)
            return false
        end
    end

    # Loads that follow the deformation (spatial tractions under finite strains) are not.
    for bc in problem.boundary_conditions
        if is_follower_load(bc, assumption)
            return false
        end
    end
    return true
end


"""
    is_follower_load(bc, assumption)

Return `true` when `bc` prescribes a quantity per unit current area under finite strains,
so that its pull-back (by Nanson's formula) depends on the deformation.
"""
is_follower_load(::AbstractBoundaryCondition, ::AbstractKinematicAssumption) = false

function is_follower_load(bc::NaturalBC, ::FiniteStrain)
    follower = quantity_configuration(bc_quantity(bc)) isa CurrentConfiguration
    return follower
end


"""
    FactorisationKey(Δt, parameters_hash, constrained_dofs_hash)

What a cached factorisation of the tangent depends on. A solver reuses its factorisation
for as long as the key is unchanged. Changing Dirichlet *values* does not change the key
(only the right-hand side, through Ferrite's `apply_rhs!`); changing *which* degrees of
freedom are constrained does.

# Fields
- `Δt::T`: the time step (the integrator's rate factors depend on it).
- `parameters_hash::UInt`: a hash of the material parameters.
- `constrained_dofs_hash::UInt`: a hash of the set of constrained degrees of freedom.
"""
struct FactorisationKey{T<:Real}
    Δt::T
    parameters_hash::UInt
    constrained_dofs_hash::UInt
end
