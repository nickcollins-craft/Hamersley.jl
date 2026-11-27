#=
Linearity.

A problem is linear when its residual is affine in the unknowns. That property is inferred
from the graph (every term's test side, flux and measure must be affine, including the
boundary terms that natural and Robin conditions become), so a user's `linear = true` is
checked once, in the formulation verification, and then trusted by the solver.

Linearity says nothing about internal state. A linear Maxwell model is affine in the
unknowns, but its old stress ``σ_n`` must still be updated and committed after every step,
so whether any law carries state (`has_internal_state`) is a separate question that the
cached-factorisation strategy also asks.

Linearity is one of two independent properties that let a solver skip work. The other is
that the tangent has not changed since it was last factorised, which is what
`FactorisationKey` records. A linear problem with a fixed time step needs one factorisation
for the whole analysis; a linear problem with adaptive time steps refactorises whenever the
step changes but never iterates; and a nonlinear problem can reuse a factorisation (modified
Newton) without being linear.
=#


"""
    is_affine(problem)
    is_affine(problem, phase)

Return `true` when the residual of `problem` is an affine function of its unknowns, in
every phase or in the one phase `phase` (a name or a position).

Boundary conditions take part through the boundary terms they become: a follower load is
integrated over the current boundary, whose measure depends on the motion under finite
strains, so it makes its phase nonlinear without any special case.
"""
is_affine(problem::Problem) = all_terms_affine(problem, problem_terms(problem))
is_affine(problem::Problem, phase) =
    all_terms_affine(problem, problem_terms(problem, phase))


"""
    all_terms_affine(problem, pairs)

Return `true` when every term of the `(equation, term)` pairs `pairs` is affine under the
kinematic assumption of `problem`.
"""
function all_terms_affine(problem::Problem, pairs)
    assumption = kinematic_assumption(problem)

    # Every term must have an affine test side (so that its variation is constant), an
    # affine flux, and a measure that does not depend on the unknowns.
    for (_, term) in pairs
        if !is_affine(term, assumption)
            return false
        end
    end
    return true
end


"""
    has_internal_state(problem)

Return `true` when some law of `problem` carries internal state at its quadrature points
(that is, its [`initial_state`](@ref) is not `nothing`), so that every step must update
and commit that state, whichever solver strategy is used.
"""
function has_internal_state(problem::Problem)
    # One stateful law is enough for the whole problem to need a state update.
    for binding in problem_constitutive(problem)
        if initial_state(binding.law) !== nothing
            return true
        end
    end
    return false
end


"""
    FactorisationKey(Δt, parameters_hash, constrained_dofs_hash)

What a cached factorisation of the tangent depends on. A solver reuses its factorisation
for as long as the key is unchanged. Changing essential *values* does not change the key
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
